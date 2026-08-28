"""Group pairing and fine alignment (GALAXY Steps 3 and 4).

Implements Sections 2.3 (Peak Group Pairing) and 2.4 (Fine Alignment
Assessment) of Deng, Zhang & Zhang (2026).
"""

import random
import warnings

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr
from tqdm import tqdm

from .util import (
    in_range_lookup,
    comp_clusters,
    get_unk_comp_clusters,
)


class AnnDataMALDI(object):
    """Orchestrator for GALAXY's pairing and fine-alignment steps.

    Parameters
    ----------
    AnnDataUnk : anndata.AnnData
        The unknown spectrum (to be aligned).
    AnnDataRef : anndata.AnnData
        The reference spectrum.
    """

    def __init__(self, AnnDataUnk, AnnDataRef):
        self.AnnDataUnk = AnnDataUnk
        self.AnnDataRef = AnnDataRef
        self.mz_valueUnk = AnnDataUnk.var
        self.mz_valueRef = AnnDataRef.var
        self.meanspectrumUnk = np.mean(AnnDataUnk.X, axis=0)
        self.meanspectrumRef = np.mean(AnnDataRef.X, axis=0)

    def ClusterPrep(self, cluster):
        """Prepare reference / unknown cluster lists from joint peak groups.

        :param cluster: the joint peak group object (jointcluster) from PeakCalling.
        """
        mz_unk = list(self.mz_valueUnk["m/z"])
        select_list_1 = in_range_lookup(mz_unk, cluster)
        mz_ref = list(self.mz_valueRef["m/z"])
        select_list_2 = in_range_lookup(mz_ref, cluster)
        select_list = [(e_i | e_j) for (e_i, e_j) in zip(select_list_1, select_list_2)]
        #
        self.ref_clusters = comp_clusters(select_list, self.mz_valueRef.index)
        self.unk_clusters = get_unk_comp_clusters(self.ref_clusters, self.mz_valueUnk.index)
        if (len(self.ref_clusters) > len(self.unk_clusters)):
            self.ref_clusters = self.ref_clusters[0:len(self.unk_clusters)]
        self.nclusters = len(self.ref_clusters)
        #
        self.midpoints_unk = [clusteri[len(clusteri) // 2] for clusteri in self.unk_clusters]

    def get_corr_peakgroup_refined(self, cluster):
        """Compute the similarity matrix across each peak group.

        Uses Pearson's correlation with a sliding-window refinement of +/- 4 grid
        points, that is four spacings of the m/z grid rather than an m/z interval of
        four, so the physical width of the search scales with the gridding increment
        (Section 2.3 of the paper).

        Rows of the returned matrix are indexed by unknown peak groups and columns by
        reference peak groups. The matrix is square but is not symmetric, because the
        unknown window is centred on the midpoint of the unknown group while the
        reference group is read from its own boundaries. Candidate pairs are
        restricted to a band of +/- 2 groups around the diagonal; entries outside the
        band are left at zero and can never be selected.

        :param cluster: the joint peak group object (jointcluster) from PeakCalling.
        """
        self.cluster = cluster
        self.ClusterPrep(cluster)
        PearsonMatrix = np.zeros((self.nclusters, self.nclusters))
        PearsonMatrixFull = np.zeros((self.nclusters * 9, self.nclusters * 9))
        for i_unk in tqdm(range(self.nclusters)):
            for i_ref in range(max(i_unk - 2, 0), min(i_unk + 3, self.nclusters)):
                midindexint = list(self.mz_valueUnk.index).index(self.midpoints_unk[i_unk])
                headindexint_ref = list(self.mz_valueRef.index).index(self.ref_clusters[i_ref][0])
                mzdiff = self.mz_valueUnk["m/z"][midindexint] - self.mz_valueRef["m/z"][
                    list(self.mz_valueRef.index).index(self.ref_clusters[i_ref][len(self.ref_clusters[i_ref]) // 2])
                ]
                Pearson_Record = -2
                Can_align = True
                for slide_unk in range(-4, 5):
                    for slide_ref in range(-4, 5):
                        indexlist_unk = slide_unk + np.arange(
                            midindexint - len(self.ref_clusters[i_ref]) // 2,
                            midindexint + len(self.ref_clusters[i_ref]) - len(self.ref_clusters[i_ref]) // 2,
                        )
                        indexlist_ref = slide_ref + np.arange(
                            headindexint_ref,
                            headindexint_ref + len(self.ref_clusters[i_ref]),
                        )
                        if indexlist_unk[0] < 0 or indexlist_ref[0] < 0:
                            # numpy indexes negatives cyclically, so the IndexError
                            # below only catches the high-m/z end of the spectrum
                            Can_align = False
                            continue
                        try:
                            specsub_unk = self.meanspectrumUnk[indexlist_unk.tolist()]
                            specsub_ref = self.meanspectrumRef[indexlist_ref.tolist()]
                            PearsonMatrixFull[i_unk * 9 + slide_unk + 4, i_ref * 9 + slide_ref + 4] = pearsonr(specsub_unk, specsub_ref)[0]
                            if PearsonMatrixFull[i_unk * 9 + slide_unk + 4, i_ref * 9 + slide_ref + 4] > Pearson_Record:
                                Pearson_Record = PearsonMatrixFull[i_unk * 9 + slide_unk + 4, i_ref * 9 + slide_ref + 4]
                        except IndexError:
                            Can_align = False
                if Can_align:
                    PearsonMatrix[i_unk, i_ref] = Pearson_Record / (abs(mzdiff) + 1)
                else:
                    PearsonMatrix[i_unk, i_ref] = -1
        self.PearsonMatrix = PearsonMatrix
        self.PearsonMatrixFull = PearsonMatrixFull

    def group_align_onestep(self, matrix, criteria, origin):
        if all(matrix.shape):
            ind = np.argwhere(matrix == matrix.max())
            maxvalue = matrix[ind[0, 0], ind[0, 1]]
            ind0 = np.array([ind[0, :]])
            if maxvalue > criteria:
                record = list(ind0 + origin)
                results1 = self.group_align_onestep(matrix[0:ind0[0, 0], 0:ind0[0, 1]], criteria, origin)
                results2 = self.group_align_onestep(matrix[(ind0[0, 0] + 1):, (ind0[0, 1] + 1):], criteria, ind0 + origin + 1)
                record = record + results1 + results2
                return record
            else:
                return []
        return []

    def peak_group_pairing(self, criteria=0):
        """GALAXY Step 3: greedily pair peak groups based on the similarity matrix.

        Operates on ``self.PearsonMatrix`` (computed by ``get_corr_peakgroup_refined``).
        Pairs survive only when the (distance-penalised) Pearson correlation
        exceeds ``criteria`` (delta in the paper). The result is stored on
        ``self.align_group`` as an n-by-2 matrix of matched (unknown, reference)
        peak-group indices.

        :param criteria: minimum (distance-penalised) Pearson correlation. Defaults to 0.
        """
        alignlist = self.group_align_onestep(self.PearsonMatrix, criteria, origin=(0, 0))
        align_results = np.empty((len(alignlist), 2))
        for i, ind in enumerate(alignlist):
            align_results[i, 0] = ind[0]
            align_results[i, 1] = ind[1]
        align_results_sort = align_results[align_results[:, 1].argsort()]
        self.align_group = align_results_sort.astype(int)
        self.align_group = np.unique(self.align_group, axis=0)
        self.nclusters = self.align_group.shape[0]

    def greedy_match(self, criteria=0):
        warnings.warn(
            "greedy_match() is deprecated; use peak_group_pairing() to match the manuscript terminology.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.peak_group_pairing(criteria=criteria)

    def fine_alignment_assessment(self, threshold=0.2, ignore=False, **kwargs):
        """GALAXY Step 4: rigid-body translation of unknown m/z values per paired group.

        For each paired peak group, find the +/- 4 grid-point offset whose
        diagonal-mean Pearson correlation is maximal; if that maximum exceeds
        ``threshold`` (or if ``ignore=True``), apply that integer-index shift. An
        offset of zero is a legitimate outcome: a group that is already aligned
        attains its maximum there and is left in place.

        :param threshold: minimum diagonal-mean Pearson correlation to accept the shift.
            Defaults to 0.2.
        :param ignore: if True, accept the best shift regardless of threshold.
        """
        # Backward-compat: accept the historical typo "threshould=" as an alias.
        if "threshould" in kwargs:
            warnings.warn(
                "fine_alignment_assessment(threshould=...) is deprecated; use threshold=... instead.",
                DeprecationWarning,
                stacklevel=2,
            )
            threshold = kwargs.pop("threshould")
        if kwargs:
            raise TypeError(f"Unexpected keyword arguments: {list(kwargs)}")

        aligned_mz_clusters_ref = []
        aligned_mz_clusters_unk = []
        changerecord = []
        #
        for i in tqdm(range(self.align_group.shape[0])):
            submatrix = self.PearsonMatrixFull[
                (self.align_group[i, 0] * 9):(self.align_group[i, 0] * 9 + 9),
                (self.align_group[i, 1] * 9):(self.align_group[i, 1] * 9 + 9),
            ]
            diagnallist = [np.nanmean(np.diag(submatrix, k=i)) for i in range(-4, 5)]
            diagnallist = [-1 if value != value else value for value in diagnallist]  # NaN -> -1
            if (max(diagnallist) > threshold) | ignore:
                change = -4 + diagnallist.index(max(diagnallist))
                midindexint = list(self.mz_valueUnk.index).index(self.midpoints_unk[self.align_group[i, 0]])
                headindexint_ref = list(self.mz_valueRef.index).index(self.ref_clusters[self.align_group[i, 1]][0])
                indexlist_ref = np.arange(
                    headindexint_ref,
                    headindexint_ref + len(self.ref_clusters[self.align_group[i, 1]]),
                )
                indexlist_unk = - change + np.arange(
                    midindexint - len(self.ref_clusters[self.align_group[i, 0]]) // 2,
                    midindexint + len(self.ref_clusters[self.align_group[i, 0]]) - len(self.ref_clusters[self.align_group[i, 0]]) // 2,
                )
                aligned_mz_clusters_ref.append(indexlist_ref)
                aligned_mz_clusters_unk.append(indexlist_unk)
                changerecord.append(- change)
        self.aligned_mz_clusters_unk = aligned_mz_clusters_unk
        self.aligned_mz_clusters_ref = aligned_mz_clusters_ref
        self.changerecord = changerecord

    def fine_align(self, threshould=0.2, ignore=False):
        warnings.warn(
            "fine_align() is deprecated; use fine_alignment_assessment() to match the manuscript terminology.",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.fine_alignment_assessment(threshold=threshould, ignore=ignore)

    def summarize(self):
        """Concatenate paired unknown/reference index arrays into flat alignment vectors."""
        self.unknownalign = np.concatenate(self.aligned_mz_clusters_unk)
        self.referenalign = np.concatenate(self.aligned_mz_clusters_ref)


class MALDI_SIM(object):
    def __init__(self, AnnDataMALDI):
        self.origindata = AnnDataMALDI.AnnDataRef.copy()
        self.aligned_mz_clusters_ref = AnnDataMALDI.aligned_mz_clusters_ref.copy()
        self.changerecord = AnnDataMALDI.changerecord.copy()
        self.shiftplot_data()

    def shiftplot_data(self):
        self.mz_valueRef = self.origindata.var.copy()
        shiftdata = np.zeros((len(self.aligned_mz_clusters_ref), 3))
        for i in tqdm(range(len(self.aligned_mz_clusters_ref))):
            shiftdata[i, 0] = self.mz_valueRef["m/z"][self.aligned_mz_clusters_ref[i]][0]
            shiftdata[i, 1] = self.mz_valueRef["m/z"][self.aligned_mz_clusters_ref[i]][-1]
            shiftdata[i, 2] = self.changerecord[i]
        shiftdatadf = pd.DataFrame(shiftdata)
        shiftdatadf = shiftdatadf.rename(columns={0: 'mzvalue_start', 1: 'mzvalue_end', 2: 'shift_fromorigin'})
        shiftdatadf["shift_fromprev"] = np.concatenate((
            [shiftdatadf["shift_fromorigin"][0]],
            np.diff(shiftdatadf["shift_fromorigin"]),
        ))
        self.shiftdatadf = shiftdatadf

    def addin(self, add_at_mz, addnumber):
        nearest_idx = np.where(abs(self.mz_valueRef - add_at_mz) == abs(self.mz_valueRef - add_at_mz).min())[0].max()
        for i in range(addnumber):
            self.arraydata = np.insert(self.arraydata, nearest_idx + 1, self.arraydata[:, nearest_idx], axis=1)
            self.mz_valueRef = np.insert(self.mz_valueRef, nearest_idx + 1, self.mz_valueRef[nearest_idx] + (i + 1) * 0.0001)

    def delout(self, del_at_mz, delnumber):
        for i in range(delnumber):
            nearest_idx = np.where(abs(self.mz_valueRef - del_at_mz) == abs(self.mz_valueRef - del_at_mz).min())[0].max()
            self.arraydata = np.delete(self.arraydata, nearest_idx, axis=1)
            self.mz_valueRef = np.delete(self.mz_valueRef, nearest_idx)

    def get_at_mz(self, idx):
        if idx > 0:
            return (self.shiftdatadf.iloc[idx - 1, 1] + self.shiftdatadf.iloc[idx, 0]) / 2
        else:
            return self.shiftdatadf.iloc[idx, 0] / 2

    def region_shuffle(self, ary, nregion):
        arylength = len(ary)
        intervals = list((np.array(range(1, nregion)) * arylength / nregion).astype("int"))
        if nregion > 1:
            listary = np.split(ary, intervals, axis=0)
        else:
            listary = [ary]
        [random.shuffle(listaryi) for listaryi in listary]
        newary = np.concatenate(listary)
        return newary

    def SIMULATEdata(self, shuffle, sigma, nregion):
        self.mz_valueRef = np.array(self.origindata.var["m/z"]).copy()
        self.arraydata = self.origindata.X.copy()
        unitdiff = self.mz_valueRef[1] - self.mz_valueRef[0]
        shift_resample = np.array(self.shiftdatadf["shift_fromprev"]).copy()
        if shuffle is True:
            shift_resample = self.region_shuffle(shift_resample, nregion=nregion)
        ## Shifting the unit
        for i in tqdm(range(self.shiftdatadf.shape[0])):
            if shift_resample[i] > 0:
                self.addin(add_at_mz=self.get_at_mz(i), addnumber=int(shift_resample[i]))
            if shift_resample[i] < 0:
                self.delout(del_at_mz=self.get_at_mz(i), delnumber=int(abs(shift_resample[i])))
        ## Adding Noises
        noise = np.random.lognormal(mean=0.0, sigma=sigma, size=self.arraydata.shape[0] * self.arraydata.shape[1])
        noise.shape = self.arraydata.shape
        self.newarray = self.arraydata + noise
        self.truemz = self.mz_valueRef
        self.newmz = np.array(range(0, len(self.truemz))) * unitdiff + self.truemz[0]
        self.currentshift = shift_resample
        self.shiftdatadf["currentshift_fromprev"] = shift_resample
        self.shiftdatadf["currentshift_fromorg"] = np.cumsum(shift_resample)

    def getAnnSim(self, shuffle=False, sigma=0.1, nregion=1):
        self.SIMULATEdata(shuffle, sigma, nregion)
        mz_value = pd.DataFrame(list(self.newmz), index=list(self.newmz)).astype('float')
        mz_value = mz_value.rename(columns={0: "m/z"})
        MALDISimData = sc.AnnData(X=self.newarray, var=mz_value, obs=self.origindata.obs)
        return MALDISimData

    def MSE(self, unkidx, refidx):
        refmz = np.array(self.origindata.var["m/z"]).copy()
        MSE = (np.square(self.truemz[unkidx] - refmz[refidx]))
        return [MSE, self.truemz[unkidx], refmz[refidx]]


class PGmzalign(MALDI_SIM):
    def __init__(self, AnnDataMALDI):
        self.unknowndata = AnnDataMALDI.AnnDataUnk.copy()
        self.aligned_mz_clusters_unk = AnnDataMALDI.aligned_mz_clusters_unk.copy()
        super(PGmzalign, self).__init__(AnnDataMALDI)

    def shiftplot_data(self):
        self.mz_valueRef = self.origindata.var.copy()
        shiftdata = np.zeros((len(self.aligned_mz_clusters_ref), 4))
        for i in tqdm(range(len(self.aligned_mz_clusters_ref))):
            shiftdata[i, 0] = self.mz_valueRef["m/z"][self.aligned_mz_clusters_ref[i]][0]
            shiftdata[i, 1] = self.mz_valueRef["m/z"][self.aligned_mz_clusters_ref[i]][-1]
            shiftdata[i, 2] = self.changerecord[i]
            shiftdata[i, 3] = self.aligned_mz_clusters_ref[i][0] - self.aligned_mz_clusters_unk[i][0]
        shiftdatadf = pd.DataFrame(shiftdata)
        shiftdatadf = shiftdatadf.rename(columns={
            0: 'mzvalue_start',
            1: 'mzvalue_end',
            2: 'shift_fromorigin',
            3: 'shift_unit_fromorigin',
        })
        shiftdatadf["shift_fromprev"] = np.concatenate((
            [shiftdatadf["shift_fromorigin"][0]],
            np.diff(shiftdatadf["shift_fromorigin"]),
        ))
        shiftdatadf["shift_unit_fromprev"] = np.concatenate((
            [shiftdatadf["shift_unit_fromorigin"][0]],
            np.diff(shiftdatadf["shift_unit_fromorigin"]),
        ))
        self.shiftdatadf = shiftdatadf

    def SIMULATEdata(self):
        self.mz_valueRef = np.array(self.unknowndata.var["m/z"]).copy()
        self.arraydata = self.unknowndata.X.copy()
        shift_resample = np.array(self.shiftdatadf["shift_unit_fromprev"]).copy()
        ## Shifting the unit
        for i in tqdm(range(self.shiftdatadf.shape[0])):
            if shift_resample[i] > 0:
                self.addin(add_at_mz=self.get_at_mz(i), addnumber=int(abs(shift_resample[i])))
            if shift_resample[i] < 0:
                self.delout(del_at_mz=self.get_at_mz(i), delnumber=int(abs(shift_resample[i])))
        self.newmz = np.array(self.origindata.var["m/z"]).copy()
        self.currentshift = shift_resample
        self.shiftdatadf["currentshift_fromprev"] = shift_resample
        self.shiftdatadf["currentshift_fromorg"] = np.cumsum(shift_resample)

    def getAnnSim(self):
        self.SIMULATEdata()
        mz_value = pd.DataFrame(list(self.newmz), index=list(self.newmz)).astype('float')
        mz_value = mz_value.rename(columns={0: "m/z"})
        dim = min(self.arraydata.shape[1], mz_value.shape[0])
        MALDISimData = sc.AnnData(
            X=self.arraydata[:, range(dim)],
            var=mz_value.iloc[range(dim), :],
            obs=self.unknowndata.obs,
        )
        return MALDISimData
