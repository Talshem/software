from Bio.Seq import Seq
from nupack import *
import pandas as pd
from joblib import Parallel, delayed
from fuzzysearch import find_near_matches
from typing import List, Tuple
import time
from ViennaRNA import RNA

## Hello World

config.threads = 8

KOZAK_SEQUENCE = "ACCAAAAUG"
RBS_SEQUENCE = "AGAGGAGA"
STOP_CODONS = ['UGA', 'UAG', 'UAA']
FORBIDDEN_PATTERNS = ["AAAA", "CCCC", "GGGG", "UUUU", "KKKKKK", "MMMMMM", "RRRRRR", "SSSSSS", "WWWWWW", "YYYYYY"]
START_CODON = "AUG"


class SwitchGenerator:

    SEED = 42
    VALID_NUCLEOTIDES = {'A', 'C', 'G', 'U'}

    def __init__(self, translated_gene_sequence: str, leader_sequence: str = "",
                 a_domain_size: int = 8, b_domain_size: int = 15, extra_stem_length: int = 0,
                 loop_size: int = 12, linker_pattern: str = "",
                 nucleotides_count_to_trim_gene_after: int = 50, is_kozak_down: bool = True):
        self.model = Model(material='rna', celsius=37)
        self.is_kozak_down = is_kozak_down
        self.translated_gene_sequence = translated_gene_sequence[3:nucleotides_count_to_trim_gene_after + 3]
        self.b_domain_size = b_domain_size
        self.a_domain_size = a_domain_size
        self.extra_stem_length = extra_stem_length
        self.loop_size = loop_size
        self.trigger_length = a_domain_size + b_domain_size
        self.leader = Domain(leader_sequence, name="leader", material="rna")
        self.extra_stem = Domain(f"N{extra_stem_length}", name="b end", material="rna")
        self.linker = Domain(linker_pattern, name="linker", material="rna")
        self.translated_gene = Domain(self.translated_gene_sequence, name="translated gene", material="rna")

        if is_kozak_down:
            self.loop = Domain(f"N{loop_size}", name="loop", material="rna")
            self.kozak = Domain(f"N{KOZAK_SEQUENCE}", name="kozak", material="rna")
            self.switch_structure = (f".{len(leader_sequence)}.{a_domain_size}({b_domain_size + extra_stem_length}"
                                     f".{loop_size}){b_domain_size + extra_stem_length}.{len(KOZAK_SEQUENCE) + 1}"
                                     f".{len(linker_pattern)}.{len(self.translated_gene_sequence)}")
        else:
            assert loop_size > len(KOZAK_SEQUENCE), "In Kozak up mode - loop size must be larger than the kozak length"
            self.loop = Domain(f"N{loop_size - len(KOZAK_SEQUENCE)}{KOZAK_SEQUENCE}", name="loop", material="rna")
            post_stem_offset = 3 - ((b_domain_size + extra_stem_length + len(linker_pattern)) % 3)
            self.post_stem = Domain(f"N{post_stem_offset}", name="post stem", material="rna")
            self.switch_structure = (f".{len(leader_sequence)}.{a_domain_size}({b_domain_size + extra_stem_length}"
                                     f".{loop_size}){b_domain_size + extra_stem_length}.{post_stem_offset}"
                                     f".{len(linker_pattern)}.{len(self.translated_gene_sequence)}")

        self.patterns_to_prevent = Pattern(FORBIDDEN_PATTERNS)
        self.options = DesignOptions(f_stop=0.02, seed=self.SEED, wobble_mutations=True)
        self.complex_name = 'activated_complex'

    @property
    def prescribed_stem_structure(self):
        return (f"({self.b_domain_size + self.extra_stem_length}"
                f".{self.loop_size}){self.b_domain_size + self.extra_stem_length}")

    def get_switch(self, trigger_sequence: str, healthy_sequence: str = None, is_watson_crick: bool = False) -> Tuple[
        str, float, float, str]:

        if len(trigger_sequence) != self.trigger_length:
            raise ValueError(f"The trigger sequence must be exactly {self.trigger_length} nucleotides long for the eukaryotic/human switch, but you provided {len(trigger_sequence)} nucleotides.")

        if not set(trigger_sequence).issubset(self.VALID_NUCLEOTIDES):
            invalid_chars = set(trigger_sequence) - self.VALID_NUCLEOTIDES
            raise ValueError(f"The trigger sequence contains invalid characters: {invalid_chars}. Only 'A', 'C', 'G', and 'U' are allowed.")

        assert len(trigger_sequence) == self.trigger_length

        trigger_domain = Domain(trigger_sequence, name="trigger domain", material='rna')
        trigger = TargetStrand([trigger_domain], name="trigger")
        trigger_structure = "." * len(trigger_sequence)
        trigger_complex = TargetComplex([trigger], trigger_structure, name='trigger complex')

        if is_watson_crick:
            trigger_reverse_complement = str(Seq(trigger_sequence).reverse_complement_rna())
            a_domain = Domain(f"{trigger_reverse_complement[:self.a_domain_size]}", name="a", material="rna")
            b_domain = Domain(f"{trigger_reverse_complement[-self.b_domain_size:]}", name="b", material="rna")
        else:
            a_domain = Domain(f"N{self.a_domain_size}", name="a", material="rna")
            b_domain = Domain(f"N{self.b_domain_size}", name="b", material="rna")

        if self.is_kozak_down:
            switch = TargetStrand([
                self.leader, a_domain, b_domain, self.extra_stem,
                self.loop, ~self.extra_stem, ~b_domain, self.kozak, self.linker, self.translated_gene
            ], name="switch")
        else:
            switch = TargetStrand([
                self.leader, a_domain, b_domain, self.extra_stem,
                self.loop, ~self.extra_stem, ~b_domain, self.post_stem, self.linker, self.translated_gene
            ], name="switch")
            stop_codon_pattern = Pattern(STOP_CODONS, scope=[~self.extra_stem, ~b_domain, self.post_stem])

        switch_complex = TargetComplex([switch], self.switch_structure, name='switch complex')

        if healthy_sequence:
            healthy_sequence_domain = Domain(healthy_sequence, name="healthy trigger domain", material='rna')
            healthy_trigger = TargetStrand([healthy_sequence_domain], name="healthy trigger")
            healthy_trigger_structure = "." * len(healthy_sequence)
            healthy_trigger_complex = TargetComplex([healthy_trigger], healthy_trigger_structure,
                                                    name='healthy trigger complex')

        activated_structure = self._get_activated_structure(switch)
        activated_complex = TargetComplex([trigger, switch], activated_structure, name=self.complex_name)

        reactants = TargetTube(on_targets={switch_complex: 1e-08, trigger_complex: 1e-08},
                               off_targets=SetSpec(max_size=2, exclude=tuple([activated_complex])), name='reactants')

        if healthy_sequence:
            off_target_reactants = TargetTube(on_targets={switch_complex: 1e-08, healthy_trigger_complex: 1e-08},
                                              off_targets=SetSpec(max_size=2), name='off_target_reactants')

        products = TargetTube(on_targets={activated_complex: 1e-08}, off_targets=SetSpec(max_size=2), name='products')

        if healthy_sequence:
            tubes = [reactants, products, off_target_reactants]
        else:
            tubes = [reactants, products]

        my_design = tube_design(tubes=tubes,
                                hard_constraints=[] if is_watson_crick or self.is_kozak_down else [stop_codon_pattern],
                                soft_constraints=[],
                                model=self.model, options=self.options)

        design_results = my_design.run(trials=1)[0]

        switch_designed_strand = [v for k, v in design_results.to_analysis.strands.items() if
                                  v.name == "switch"].pop()
        ensemble_defect = self._get_ensemble_defect(design_results)
        complex_concentration = self._get_trigger_switch_complex_concentration(design_results, self.complex_name)
        switch_mfe_structure = self._get_switch_mfe_structure(design_results, switch_complex)
        return str(switch_designed_strand), ensemble_defect, complex_concentration, switch_mfe_structure

    def _get_activated_structure(self, switch) -> str:
        return ('(' * self.trigger_length + '+' + '.' * self.leader.nt() + ')' * self.trigger_length +
                '.' * (switch.nt() - self.leader.nt() - self.trigger_length))

    def _get_ensemble_defect(self, design_result) -> float:
        return design_result.defects.ensemble_defect

    def _get_trigger_switch_complex_concentration(self, design_result, complex_name) -> float:
        complex_concentrations = design_result.concentrations.table[
            design_result.concentrations.table.complex_name == complex_name].iloc[0]
        target_concentration = complex_concentrations["target_concentration"]
        predicted_concentration = complex_concentrations["concentration"]
        return predicted_concentration / target_concentration

    def _get_switch_mfe_structure(self, design_result, switch_target_strand) -> str:
        return str(mfe(strands=design_result.to_analysis(switch_target_strand), model=self.model)[0].structure)


class ProkaryoticSwitchGenerator:

    SEED = 42
    VALID_NUCLEOTIDES = {'A', 'C', 'G', 'U'}

    def __init__(self, translated_gene_sequence: str, leader_sequence: str = "GGG",
                 a_domain_size: int = 12, b_domain_size_pre_bulge: int = 9, b_domain_size_post_bulge: int = 6,
                 loop_size: int = 11, linker_pattern: str = "AACCUGGCGGCAGCGCAAAAG",
                 nucleotides_count_to_trim_gene_after: int = 30):
        self.model = Model(material='rna', celsius=37)
        self.translated_gene_sequence = translated_gene_sequence[3:nucleotides_count_to_trim_gene_after + 3]
        self.b_domain_size_pre_bulge = b_domain_size_pre_bulge
        self.b_domain_size_post_bulge = b_domain_size_post_bulge
        self.a_domain_size = a_domain_size
        self.loop_size = loop_size
        self.trigger_length = a_domain_size + b_domain_size_pre_bulge + b_domain_size_post_bulge + len(START_CODON)
        self.leader = Domain(leader_sequence, name="leader", material="rna")
        self.linker = Domain(linker_pattern, name="linker", material="rna")
        self.translated_gene = Domain(self.translated_gene_sequence, name="translated gene", material="rna")

        assert loop_size > len(RBS_SEQUENCE), "loop size must be larger than the RBS length"

        self.loop = Domain(f"N{loop_size - len(RBS_SEQUENCE)}{RBS_SEQUENCE}", name="loop", material="rna")
        self.bulge = Domain(f"{START_CODON}", name="bulge", material="rna")
        self.switch_structure = (f".{len(leader_sequence)}.{a_domain_size}({b_domain_size_pre_bulge}"
                                 f".{len(START_CODON)}({b_domain_size_post_bulge}"
                                 f".{loop_size}){b_domain_size_post_bulge}.{len(START_CODON)}"
                                 f"){b_domain_size_pre_bulge}"
                                 f".{len(linker_pattern)}.{len(self.translated_gene_sequence)}")

        self.patterns_to_prevent = Pattern(FORBIDDEN_PATTERNS)
        self.options = DesignOptions(f_stop=0.02, seed=self.SEED, wobble_mutations=True)
        self.complex_name = 'activated_complex'

    @property
    def prescribed_stem_structure(self):
        return (f"({self.b_domain_size_pre_bulge}"
                f".{len(START_CODON)}({self.b_domain_size_post_bulge}"
                f".{self.loop_size}){self.b_domain_size_post_bulge}.{len(START_CODON)}"
                f"){self.b_domain_size_pre_bulge}")

    def get_switch(self, trigger_sequence: str, healthy_sequence: str = None,
                   is_watson_crick: bool = True) -> Tuple[str, float, float, str]:

        if len(trigger_sequence) != self.trigger_length:
            raise ValueError(f"The trigger sequence must be exactly {self.trigger_length} nucleotides long for the prokaryotic switch, but you provided {len(trigger_sequence)} nucleotides.")

        if not set(trigger_sequence).issubset(self.VALID_NUCLEOTIDES):
            invalid_chars = set(trigger_sequence) - self.VALID_NUCLEOTIDES
            raise ValueError(f"The trigger sequence contains invalid characters: {invalid_chars}. Only 'A', 'C', 'G', and 'U' are allowed.")

        assert len(trigger_sequence) == self.trigger_length

        trigger_domain = Domain(trigger_sequence, name="trigger domain", material='rna')
        trigger = TargetStrand([trigger_domain], name="trigger")
        trigger_structure = "." * len(trigger_sequence)
        trigger_complex = TargetComplex([trigger], trigger_structure, name='trigger complex')

        if is_watson_crick:
            trigger_reverse_complement = str(Seq(trigger_sequence).reverse_complement_rna())
            a_domain = Domain(f"{trigger_reverse_complement[:self.a_domain_size]}", name="a", material="rna")
            b_domain_pre_bulge = Domain(
                f"{trigger_reverse_complement[self.a_domain_size:self.a_domain_size + self.b_domain_size_pre_bulge]}",
                name="b_domain_pre_bulge", material="rna")
            b_domain_bulge = Domain(
                f"{trigger_reverse_complement[self.a_domain_size + self.b_domain_size_pre_bulge:self.a_domain_size + self.b_domain_size_pre_bulge + len(START_CODON)]}",
                name="b_domain_bulge", material="rna")
            b_domain_post_bulge = Domain(
                f"{trigger_reverse_complement[-self.b_domain_size_post_bulge:]}",
                name="b_domain_post_bulge", material="rna")
        else:
            a_domain = Domain(f"N{self.a_domain_size}", name="a", material="rna")
            b_domain_pre_bulge = Domain(f"N{self.b_domain_size_pre_bulge}",
                                        name="b_domain_pre_bulge", material="rna")
            b_domain_bulge = Domain(f"N{len(START_CODON)}", name="b_domain_bulge", material="rna")
            b_domain_post_bulge = Domain(f"N{self.b_domain_size_post_bulge}",
                                         name="b_domain_post_bulge", material="rna")

        switch = TargetStrand([
            self.leader, a_domain, b_domain_pre_bulge, b_domain_bulge, b_domain_post_bulge,
            self.loop, ~b_domain_post_bulge, self.bulge, ~b_domain_pre_bulge, self.linker, self.translated_gene
        ], name="switch")
        stop_codon_pattern = Pattern(STOP_CODONS, scope=[~b_domain_pre_bulge])

        switch_complex = TargetComplex([switch], self.switch_structure, name='switch complex')

        if healthy_sequence:
            healthy_sequence_domain = Domain(healthy_sequence, name="healthy trigger domain", material='rna')
            healthy_trigger = TargetStrand([healthy_sequence_domain], name="healthy trigger")
            healthy_trigger_structure = "." * len(healthy_sequence)
            healthy_trigger_complex = TargetComplex([healthy_trigger], healthy_trigger_structure,
                                                    name='healthy trigger complex')

        activated_structure = self._get_activated_structure(switch)
        activated_complex = TargetComplex([trigger, switch], activated_structure, name=self.complex_name)

        reactants = TargetTube(on_targets={switch_complex: 1e-08, trigger_complex: 1e-08},
                               off_targets=SetSpec(max_size=2, exclude=tuple([activated_complex])), name='reactants')

        if healthy_sequence:
            off_target_reactants = TargetTube(on_targets={switch_complex: 1e-08, healthy_trigger_complex: 1e-08},
                                              off_targets=SetSpec(max_size=2), name='off_target_reactants')

        products = TargetTube(on_targets={activated_complex: 1e-08}, off_targets=SetSpec(max_size=2), name='products')

        if healthy_sequence:
            tubes = [reactants, products, off_target_reactants]
        else:
            tubes = [reactants, products]

        my_design = tube_design(tubes=tubes,
                                hard_constraints=[],
                                soft_constraints=[],
                                model=self.model, options=self.options)

        design_results = my_design.run(trials=1)[0]

        switch_designed_strand = [v for k, v in design_results.to_analysis.strands.items() if
                                  v.name == "switch"].pop()
        ensemble_defect = self._get_ensemble_defect(design_results)
        complex_concentration = self._get_trigger_switch_complex_concentration(design_results, self.complex_name)
        switch_mfe_structure = self._get_switch_mfe_structure(design_results, switch_complex)
        return str(switch_designed_strand), ensemble_defect, complex_concentration, switch_mfe_structure

    def _get_activated_structure(self, switch) -> str:
        return ('(' * self.trigger_length + '+' + '.' * self.leader.nt() + ')' * self.trigger_length +
                '.' * (switch.nt() - self.leader.nt() - self.trigger_length))

    def _get_ensemble_defect(self, design_result) -> float:
        return design_result.defects.ensemble_defect

    def _get_trigger_switch_complex_concentration(self, design_result, complex_name) -> float:
        complex_concentrations = design_result.concentrations.table[
            design_result.concentrations.table.complex_name == complex_name].iloc[0]
        target_concentration = complex_concentrations["target_concentration"]
        predicted_concentration = complex_concentrations["concentration"]
        return predicted_concentration / target_concentration

    def _get_switch_mfe_structure(self, design_result, switch_target_strand) -> str:
        return str(mfe(strands=design_result.to_analysis(switch_target_strand), model=self.model)[0].structure)


# Utility function
def check_forbidden_patterns(sequence):
    """Checks if a sequence contains any forbidden patterns."""
    for pattern in FORBIDDEN_PATTERNS:
        if pattern in sequence:
            return True
    return False


class HomologySwitchGenerator:
    def __init__(self, cell_type: str, reporter_gene: str, transcripts_dict: dict):

        self.cell_type = cell_type
        self.reporter_gene = reporter_gene
        self.transcripts_dict = transcripts_dict

    def generate_switch_with_homology(self, triggers: List[str], transcripts_dict: dict,
                                      n_jobs: int = 4) -> pd.DataFrame:
        """
        Generates toehold switches with homologous sequence ranking.

        :param triggers: List of trigger sequences
        :param transcripts_dict: Genome data for homology search
        :param n_jobs: Number of parallel jobs to run
        :return: DataFrame with ranked homologous sequences and switch results
        """
        # Homology search with parallelization
        s_h = time.time()
        homo_res = Parallel(n_jobs=n_jobs)(
            delayed(self.find_homology)(trigger, transcripts_dict) for trigger in triggers)
        e_h = time.time()
        print(f'Homology search time= {e_h - s_h} seconds, n_triggers = {len(triggers)}, cell= {self.cell_type}')

        # Extract top homology sequences with parallelization
        rrf_ranks = self.extract_top_homology_sequences(homo_res)
        homology_sequences = [ranked_df['sequence'].get(0) for ranked_df in rrf_ranks]

        # Generate switches with parallelization
        s_switch = time.time()
        switch_res = Parallel(n_jobs=n_jobs)(delayed(self.generate_switch)(trigger, top_homology_sequence)
                                             for trigger, top_homology_sequence in zip(triggers, homology_sequences)
                                             )
        e_switch = time.time()
        print(f'Switch generation time= {e_switch - s_switch} seconds')

        # Combine results
        results = pd.DataFrame(switch_res, columns=['switch', 'complex_concentration'])
        results['trigger_window'] = triggers
        return results

    def find_homology(self, trigger: str, genome_data: dict) -> List[dict]:
        """
        Searches for homologous sequences in the genome data.
        """
        genes_sub_sequences = []
        for gene_data_dict in genome_data:
            gene_sub_seqs = {}
            mRNA_seq = gene_data_dict['sequence']
            gene_name = gene_data_dict['gene']
            protein = gene_data_dict['protein']
            seq_match_mapping = self.build_homology_map(trigger, mRNA_seq, gene_name, protein)
            if seq_match_mapping:
                gene_sub_seqs[gene_name] = seq_match_mapping
                genes_sub_sequences.append(gene_sub_seqs)
        return genes_sub_sequences

    def build_homology_map(self, trigger: str, seq: str, gene_name: str, protein_name: str) -> List[dict]:
        """
        Builds a homology map of the trigger against a gene sequence.
        """
        sequence_match_mapping = []
        matches = find_near_matches(trigger, seq, max_insertions=0, max_deletions=0, max_l_dist=4)
        for match_obj in matches:
            locus_start = match_obj.start
            locus_end = match_obj.end
            sub_seq = seq[locus_start:locus_end]
            seq_location = (locus_start, locus_end)
            distance = match_obj.dist
            sequence_match_mapping.append({
                'distance': distance, 'idx': seq_location,
                'sequence': sub_seq, 'gene': gene_name, 'protein': protein_name
            })
        return sequence_match_mapping

    def extract_top_homology_sequences(self, triggers_homology_mapping: List[List[dict]]) -> List[pd.DataFrame]:
        """
        Extracts the top-ranked homologous sequences for each trigger.
        """
        homo_dfs = []
        for trigger_homology in triggers_homology_mapping:
            if trigger_homology:
                trigger_all_matches = []
                for gene_homology in trigger_homology:
                    trigger_all_matches.extend(gene_homology.values())
                matches_df = pd.DataFrame(trigger_all_matches[0])

                mfe_dict = {'homologous_trigger_mfe': []}
                for homos in trigger_all_matches[0]:
                    homo_trigger = homos['sequence']
                    structure, mfe = RNA.fold(homo_trigger)
                    mfe_dict['homologous_trigger_mfe'].append(mfe)
                mfe_df = pd.DataFrame(mfe_dict)

                trig_res_df = pd.concat([matches_df, mfe_df], axis=1)
            else:
                trig_res_df = pd.DataFrame(
                    columns=['distance', 'idx', 'sequence', 'gene', 'protein', 'homologous_trigger_mfe'])
            homo_dfs.append(trig_res_df)

        higher = ['homologous_trigger_mfe']
        lower = ['distance']
        rrf_rank = [self.RRF(trigger_homology_df, higher, lower, index='sequence') for trigger_homology_df in homo_dfs]
        return rrf_rank

    def RRF(self, ranking_df: pd.DataFrame, higher_is_better_cols: List[str], lower_is_better_cols: List[str],
            index: str, k: int = 60) -> pd.DataFrame:
        """
        Rank fusion (RRF) algorithm to rank sequences based on multiple metrics.
        """
        ranking_df = ranking_df.copy().reset_index(drop=True)

        for col in higher_is_better_cols:
            ranking_df[col + '_rank'] = ranking_df[col].rank(ascending=False)
        for col in lower_is_better_cols:
            ranking_df[col + '_rank'] = ranking_df[col].rank(ascending=True)

        ranked_columns = higher_is_better_cols + lower_is_better_cols
        ranked_columns = [f'{col}_rank' for col in ranked_columns]
        ranking_df[f'{index}_RRF'] = ranking_df[ranked_columns].apply(lambda row: sum(1 / (k + rank) for rank in row),
                                                                      axis=1)

        return ranking_df.sort_values(by=f'{index}_RRF', ascending=True)

    def generate_switch(self, trigger: str, homologous_sequence: str) -> Tuple[str, float]:
        """
        Generates a switch for the given trigger and homologous sequence.
        """
        if len(trigger) == 30:  # Prokaryotic trigger
            switch_generator = ProkaryoticSwitchGenerator(translated_gene_sequence=self.reporter_gene)
        elif len(trigger) == 23:  # Eukaryotic trigger
            switch_generator = SwitchGenerator(translated_gene_sequence=self.reporter_gene)
        else:
            raise ValueError(f"Invalid trigger length {len(trigger)}. It should be 23 for eukaryotic or 30 for prokaryotic.")

        switch_designed_strand, ensemble_defect, complex_concentration, switch_mfe_structure = switch_generator.get_switch(
            trigger_sequence=trigger, healthy_sequence=homologous_sequence)
        return switch_designed_strand, complex_concentration


if __name__ == "__main__":
    # Example genome data to simulate transcript sequences (you can replace this with real data)
    transcripts_dict = [
        {
            'gene': 'GeneA',
            'protein': 'ProteinA',
            'sequence': 'AUGGCCUAGCGCUAUGCCCUAUGGGAUGCUUCGGAUAG'
        },
        {
            'gene': 'GeneB',
            'protein': 'ProteinB',
            'sequence': 'AUGGCUCGAUUGCCCUUCUAGGAUCGUAGCUAGGAUC'
        }
    ]

    # Example list of triggers
    # First one is for eukaryotic (23 nucleotides), second one is for prokaryotic (30 nucleotides)
    triggers = [
        "AUGGCCUAGCGCUAUGCCCUAGG",  # Eukaryotic trigger (23 nucleotides)
        "AUGGCUCGAUUGCCCUUCUAGGAUCGUAGC"  # Prokaryotic trigger (30 nucleotides)
    ]

    # Example reporter gene sequence (this is the sequence to which the switch is tied)
    reporter_gene = "AUGGCCUAGCGCUAUGCCCUAUGGGAUGCUUCGGAUAG"

    # Initialize the HomologySwitchGenerator class for Homo sapiens
    homology_switch_generator = HomologySwitchGenerator(
        cell_type="Homo sapiens",  # Eukaryotic cell type
        reporter_gene=reporter_gene,  # The reporter gene sequence
        transcripts_dict=transcripts_dict  # Mock genome data for homology search
    )

    # Generate switches with homology search and ranking
    results = homology_switch_generator.generate_switch_with_homology(triggers,
                                                                      homology_switch_generator.transcripts_dict)

    # Display the results
    print(results)
