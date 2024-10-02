import re
import subprocess
import sys
from typing import List

import lightgbm
import numpy as np
import pandas as pd
from Bio.SeqUtils import gc_fraction
from RNA import RNA
from nupack import *

#sys.path.insert(1, '/Users/netanelerlich/Downloads/PyFeat-master/Codes')
sys.path.insert(1, '/workspace/PyFeat-master/Codes')
import generateFeatures
config.threads = 4

#RNACOFOLD_BINARY_NAME = '/Users/netanelerlich/ViennaRNA-2.6.4/src/bin/RNAcofold'
#RNAUP_BINARY_NAME = "/Users/netanelerlich/ViennaRNA-2.6.4/src/bin/RNAup"

RNAUP_BINARY_NAME = "/workspace/miniconda3/envs/myenv/bin/RNAup"
RNACOFOLD_BINARY_NAME = '/workspace/miniconda3/envs/myenv/bin/RNAcofold'

class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

"""
prokaryotes:

ideal_stem_structure = '..............((((((((((...((((((...........))))))...))))))))))'
trigger_binding_site_start = 3
trigger_binding_site_end = 32
stem_start = 14
stem_end = 62
loop_start = 33
loop_end = 43
stem_top_start = 27
stem_top_end = 49

eukaryotes (kozak down):

ideal_stem_structure = '........(((((((((((((((............)))))))))))...'
trigger_binding_site_start = 0
trigger_binding_site_end = 22
stem_start = 8
stem_end = 49
loop_start = 23
loop_end = 34
stem_top_start = 17
stem_top_end = 40
"""
class EukaryoticScoreCalculator:
    def __init__(self, model_file, selected_features: List[str], ideal_stem_structure: str,
                 trigger_binding_site_start: int, trigger_binding_site_end: int, stem_start: int, stem_end: int,
                 loop_start: int, loop_end: int, stem_top_start: int, stem_top_end: int):
        self.ideal_stem_structure = ideal_stem_structure
        self.trigger_binding_site_start = trigger_binding_site_start
        self.trigger_binding_site_end = trigger_binding_site_end
        self.stem_start = stem_start
        self.stem_end = stem_end
        self.loop_start = loop_start
        self.loop_end = loop_end
        self.stem_top_start = stem_top_start
        self.stem_top_end = stem_top_end
        self.sample_size = 10000
        self.pyfeat_args = {
            'sequenceType': 'RNA',
            'kGap': 2,  # default is 5
            'kTuple': 3,
            'pseudoKNC': 1,
            'zCurve': 1,
            'cumulativeSkew': 1,
            'gcContent': 0,
            'atgcRatio': 1,
            'monoMono': 1,
            'monoDi': 1,
            'monoTri': 1,
            'diMono': 1,
            'diDi': 1,
            'diTri': 1,
            'triMono': 1,
            'triDi': 1
        }
        self.selected_features = selected_features
        self.model = lightgbm.Booster(model_file=model_file)

    def get_score(self, switch, trigger) -> float:

        features_df = pd.DataFrame()

        #RNAfold
        features_df[[
            "switch_mfe_structure", "switch_mfe", "switch_mfe_matches_count", "switch_mfe_ensemble_defect",
            "switch_mfe_probability", "switch_centroid_structure", "switch_centroid_energy",
            "switch_centroid_matches_count", "switch_centroid_ensemble_defect", "switch_centroid_probability",
            "switch_partition_function",
            "estimated_probability_for_ideal_stem", "switch_estimated_energy_mean", "switch_estimated_energy_std"
        ]] = [self._get_energy_properties(switch, self.ideal_stem_structure)]

        features_df[[
            "binding_site_mfe_structure", "binding_site_mfe", "binding_site_mfe_matches_count",
            "binding_site_mfe_ensemble_defect", "binding_site_mfe_probability",
            "binding_site_centroid_structure", "binding_site_centroid_energy", "binding_site_centroid_matches_count",
            "binding_site_centroid_ensemble_defect", "binding_site_centroid_probability",
            "binding_site_partition_function",
            "binding_site_estimated_energy_mean", "binding_site_estimated_energy_std"
        ]] = [self._get_energy_properties(
            switch[int(self.trigger_binding_site_start):int(self.trigger_binding_site_end) + 1])]

        features_df[[
            "stem_mfe_structure", "stem_mfe", "stem_mfe_matches_count",
            "stem_mfe_ensemble_defect", "stem_mfe_probability",
            "stem_centroid_structure", "stem_centroid_energy", "stem_centroid_matches_count",
            "stem_centroid_ensemble_defect", "stem_centroid_probability", "stem_partition_function",
            "stem_estimated_energy_mean", "stem_estimated_energy_std"
        ]] = [self._get_energy_properties(switch[int(self.stem_start):int(self.stem_end) + 1])]

        features_df[[
            "loop_mfe_structure", "loop_mfe", "loop_mfe_matches_count",
            "loop_mfe_ensemble_defect", "loop_mfe_probability",
            "loop_centroid_structure", "loop_centroid_energy", "loop_centroid_matches_count",
            "loop_centroid_ensemble_defect", "loop_centroid_probability", "loop_partition_function",
            "loop_estimated_energy_mean", "loop_estimated_energy_std"
        ]] = [self._get_energy_properties(switch[int(self.loop_start):int(self.loop_end) + 1])]

        features_df[[
            "loop_to_end_mfe_structure", "loop_to_end_mfe", "loop_to_end_mfe_matches_count",
            "loop_to_end_mfe_ensemble_defect", "loop_to_end_mfe_probability",
            "loop_to_end_centroid_structure", "loop_to_end_centroid_energy", "loop_to_end_centroid_matches_count",
            "loop_to_end_centroid_ensemble_defect", "loop_to_end_centroid_probability",
            "loop_to_end_partition_function", "loop_to_end_estimated_energy_mean",
            "loop_to_end_estimated_energy_std"
        ]] = [self._get_energy_properties(switch[int(self.loop_start):int(self.stem_end) + 1])]

        features_df[[
            "loop_to_stem_end_mfe_structure", "loop_to_stem_end_mfe", "loop_to_stem_end_mfe_matches_count",
            "loop_to_stem_end_mfe_ensemble_defect", "loop_to_stem_end_mfe_probability",
            "loop_to_stem_end_centroid_structure", "loop_to_stem_end_centroid_energy",
            "loop_to_stem_end_centroid_matches_count",
            "loop_to_stem_end_centroid_ensemble_defect", "loop_to_stem_end_centroid_probability",
            "loop_to_stem_end_partition_function", "loop_to_stem_end_estimated_energy_mean",
            "loop_to_stem_end_estimated_energy_std"
        ]] = [self._get_energy_properties(switch[int(self.loop_start):])]

        features_df[[
            "stem_start_loop_end_mfe_structure", "stem_start_loop_end_mfe", "stem_start_loop_end_mfe_matches_count",
            "stem_start_loop_end_mfe_ensemble_defect", "stem_start_loop_end_mfe_probability",
            "stem_start_loop_end_centroid_structure", "stem_start_loop_end_centroid_energy",
            "stem_start_loop_end_centroid_matches_count",
            "stem_start_loop_end_centroid_ensemble_defect", "stem_start_loop_end_centroid_probability",
            "stem_start_loop_end_partition_function", "stem_start_loop_end_estimated_energy_mean",
            "stem_start_loop_end_estimated_energy_std"
        ]] = [self._get_energy_properties(switch[int(self.stem_start):int(self.loop_end) + 1])]

        features_df[[
            "start_loop_end_mfe_structure", "start_loop_end_mfe", "start_loop_end_mfe_matches_count",
            "start_loop_end_mfe_ensemble_defect", "start_loop_end_mfe_probability",
            "start_loop_end_centroid_structure", "start_loop_end_centroid_energy",
            "start_loop_end_centroid_matches_count",
            "start_loop_end_centroid_ensemble_defect", "start_loop_end_centroid_probability",
            "start_loop_end_partition_function", "start_loop_end_estimated_energy_mean",
            "start_loop_end_estimated_energy_std"
        ]] = [self._get_energy_properties(switch[:int(self.loop_end) + 1])]

        features_df[[
            "stem_top_mfe_structure", "stem_top_mfe", "stem_top_mfe_matches_count",
            "stem_top_mfe_ensemble_defect", "stem_top_mfe_probability",
            "stem_top_centroid_structure", "stem_top_centroid_energy",
            "stem_top_centroid_matches_count",
            "stem_top_centroid_ensemble_defect", "stem_top_centroid_probability",
            "stem_top_partition_function", "stem_top_estimated_energy_mean",
            "stem_top_estimated_energy_std"
        ]] = [self._get_energy_properties(switch[int(self.stem_top_start):int(self.stem_top_end) + 1])]

        features_df[[
            "trigger_mfe_structure", "trigger_mfe", "trigger_mfe_matches_count",
            "trigger_mfe_ensemble_defect", "trigger_mfe_probability",
            "trigger_centroid_structure", "trigger_centroid_energy", "trigger_centroid_matches_count",
            "trigger_centroid_ensemble_defect", "trigger_centroid_probability", "trigger_partition_function",
            "trigger_estimated_energy_mean", "trigger_estimated_energy_std"
        ]] = [self._get_energy_properties(trigger)]

        # RNAup
        features_df[["trigger_with_switch_total_free_energy", "trigger_with_switch_duplex_formation_energy",
                     "switch_opening_energy", "triger_opening_energy", "matches_count"]] = (
            [self._get_trigger_with_switch_mfe(trigger,switch)])

        features_df[["trigger_with_trigger_total_free_energy", "trigger_with_trigger_duplex_formation_energy",
                     "trigger_opening_energy_", "trigger_opening_energy_",
                     "trigger_with_trigger_matches_count"]] = [self._get_trigger_with_switch_mfe(trigger, trigger)]

        features_df[["switch_with_switch_total_free_energy", "switch_with_switch_duplex_formation_energy",
                     "switch_opening_energy_", "trigger_opening_energy_",
                     "switch_with_switch_matches_count"]] = [self._get_trigger_with_switch_mfe(trigger, trigger)]

        features_df = features_df.drop(columns=['switch_opening_energy_', 'trigger_opening_energy_'])

        # RNAcofold
        features_df[["trigger_with_switch_mfe_cofold",
                     "trigger_with_switch_matches_count_cofold"]] = [self._get_trigger_with_switch_cofold_mfe(trigger, switch)]

        # gc_content
        features_df["switch_gc_content"] = [gc_fraction(switch)]
        features_df["binding_site_gc_content"] = [gc_fraction(
            switch[int(self.trigger_binding_site_start):int(self.trigger_binding_site_end) + 1])]
        features_df["stem_gc_content"] = [gc_fraction(switch[int(self.stem_start):int(self.stem_end) + 1])]
        features_df["loop_gc_content"] = [gc_fraction(switch[int(self.loop_start):int(self.loop_end) + 1])]
        features_df["loop_to_stem_end_gc_content"] = [gc_fraction(switch[int(self.loop_start):int(self.stem_end) + 1])]
        features_df["loop_to_end_gc_content"] = [gc_fraction(switch[int(self.loop_start):])]
        features_df["stem_start_to_loop_end_gc_content"] = [gc_fraction(
            switch[int(self.stem_start):int(self.loop_end) + 1])]
        features_df["start_to_loop_end_gc_content"] = [gc_fraction(switch[:int(self.loop_end) + 1])]
        features_df["stem_top_gc_content"] = [gc_fraction(
            switch[int(self.stem_top_start):int(self.stem_top_end) + 1])]
        features_df["trigger_gc_content"] = [gc_fraction(trigger)]

        # nupack
        features_df["nupack_predicted_concentration"] = [self._get_nupack_properties(trigger, switch)]

        features_df["is_au_before_hairpin"] = [1 if switch[self.stem_start - 1] in "AU" else 0]

        X = features_df.drop(columns=[feature for feature in features_df.columns if "_structure" in feature])
        T = generateFeatures.gF(dotdict(self.pyfeat_args), [switch], [0])
        extended_X = pd.DataFrame(np.column_stack([X, T[:, :-1]]))  # last column is the label
        pyfeat_features = [f'zCurve_{i}' for i in range(1, 4)] + [f'cumulativeSkew_{i}' for i in range(1, 3)] + [
            'atgcRatio'] + [f'pseudoKNC_{i}' for i in range(1, 1 + 4 + 16 + 64)] + [
            f'monoMonoKGap_{i}' for i in range(1, 33)] + [f'monoDiKGap_{i}' for i in range(1, 129)] + [
            f'monoTriKGap_{i}' for i in range(1, 513)] + [f'diMonoKGap_{i}' for i in range(1, 129)] + [
            f'diDiKGap_{i}' for i in range(1, 513)] + [f'diTriKGap_{i}' for i in range(1, 2049)] + [
            f'triMonoKGap_{i}' for i in range(1, 513)] + [f'triDiKGap_{i}' for i in range(1, 2049)]
        extended_X.columns = list(X.columns) + pyfeat_features
        prediction = self.model.predict(extended_X[self.selected_features].astype(np.float64))[0]
        return prediction



    def _get_energy_properties(self, sequence: str, ideal_stem_structure: str = None):
        # create a new model details structure
        md = RNA.md()

        # change temperature and dangle model
        md.noLP = 1  # 20 Deg Celcius
        md.dangles = 2  # Dangle Model 1
        md.uniq_ML = 1

        fc = RNA.fold_compound(sequence, md)
        (mfe_structure, mfe) = fc.mfe()
        fc.exp_params_rescale(mfe)
        (pp, pf) = fc.pf()
        (centroid_structure, distance) = fc.centroid()
        centroid_energy = fc.eval_structure(centroid_structure)
        centroid_ensemble_defect = fc.ensemble_defect(centroid_structure)
        mfe_ensemble_defect = fc.ensemble_defect(mfe_structure)
        centroid_probability = fc.pr_structure(centroid_structure)
        mfe_probability = fc.pr_structure(mfe_structure)
        mfe_match_count = mfe_structure.count("(")
        centroid_match_count = centroid_structure.count("(")

        sampled_structures = fc.pbacktrack(self.sample_size)
        energies = np.array([fc.eval_structure(sampled_structure) for sampled_structure in sampled_structures])

        if ideal_stem_structure:
            estimated_probability_for_ideal_stem = sum(
                [ideal_stem_structure == sampled_structure[:len(ideal_stem_structure)]
                 for sampled_structure in sampled_structures]) / self.sample_size

            return (mfe_structure, mfe, mfe_match_count, mfe_ensemble_defect, mfe_probability,
                    centroid_structure, centroid_energy, centroid_match_count, centroid_ensemble_defect,
                    centroid_probability,
                    pf, estimated_probability_for_ideal_stem, energies.mean(), energies.std())

        return (mfe_structure, mfe, mfe_match_count, mfe_ensemble_defect, mfe_probability,
                centroid_structure, centroid_energy, centroid_match_count, centroid_ensemble_defect,
                centroid_probability,
                pf, energies.mean(), energies.std())


    def _get_trigger_with_switch_mfe(self, trigger, switch) -> tuple[float, float, float, float, int]:
        rna_up_output = self._run_rna_up(trigger, switch)

        try:
            regex_result = re.search("\([0-9-.]+ = [0-9-.]+ \+ [0-9-+.]+ \+ ([0-9-+.]+\)|inf)",
                                     rna_up_output).group()
        except Exception as e:
            print(trigger)
            print(switch)
            raise e

        string_values = regex_result.replace("(", "").replace(")", "").replace("=", "").replace("+", "").split()
        values = [float(value) for value in string_values]
        total_free_energy = values[0]
        energy_from_duplex_formation = values[1]
        switch_opening_energy = values[2]
        trigger_opening_energy = values[3]
        secondary_structure = rna_up_output.split()[0]
        match_count = secondary_structure.count("(")
        return total_free_energy, energy_from_duplex_formation, switch_opening_energy, trigger_opening_energy, match_count

    def _run_rna_up(self, sequence1, sequence2):
        try:
            return subprocess.check_output(
                [RNAUP_BINARY_NAME, "-b", "-d2", "--noLP", "-c", "'S'", "-o"],
                universal_newlines=True,
                input=f"{sequence1}&{sequence2}",
                text=True
            )
        except Exception as e:
            print(f"Error while running RNAup on {sequence1} and {sequence2}. Error: {e}")
            raise

    def _get_trigger_with_switch_cofold_mfe(self, trigger, switch) -> tuple[float, int]:
        rna_cofold_output = self._run_rna_cofold(trigger, switch)
        regex_result = re.search("\([0-9-]+[0-9.]*\)", rna_cofold_output).group()
        mfe = float(regex_result.replace('(', '').replace(')', ''))
        secondary_structure = rna_cofold_output.split()[1]
        match_count = secondary_structure.count("(")
        return mfe, match_count

    def _run_rna_cofold(self, sequence1, sequence2):
        try:
            return subprocess.check_output(
                [RNACOFOLD_BINARY_NAME, "-d2", "--noLP"],
                universal_newlines=True,
                input=f"{sequence1}&{sequence2}",
                text=True
            )
        except Exception as e:
            print(f"Error while running RNAcofold on {sequence1} and {sequence2}. Error: {e}")
            raise

    def _get_nupack_properties(self, switch, trigger_strand):
        # specify strands
        trigger_strand = Strand(trigger_strand, name='trigger', material="rna")
        switch_strand = Strand(switch, name='switch', material="rna")

        # specify tubes
        tube = Tube(strands={trigger_strand: 1e-8, switch_strand: 1e-8}, complexes=SetSpec(max_size=2), name='t1')
        tube_results = tube_analysis(tubes=[tube], model=Model(material='rna', celsius=37))
        complex_concentration = self._get_trigger_switch_complex_concentration(tube_results)

        return complex_concentration

    def _get_trigger_switch_complex_concentration(self, design_result) -> float:
        target_concentration = sum([c for name, c in design_result["t1"].complex_concentrations.items()])

        try:
            predicted_concentration = [(name, c) for name, c in design_result["t1"].complex_concentrations.items()
                                       if str(name) == "<Complex (switch+trigger)>"][0][1]
        except IndexError as e:
            predicted_concentration = [(name, c) for name, c in design_result["t1"].complex_concentrations.items()
                                       if str(name) == "<Complex (trigger+switch)>"][0][1]

        return predicted_concentration / target_concentration
