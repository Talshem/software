from typing import Tuple
from viennaRNA import mfe

from nupack import *
config.threads = 8

KOZAK_SEQUENCE = {'Homo_sapiens': 'GCCACCATGGCC', 'Escherichia_coli': 'AGGAGGGGCAAATATG', 'Saccharomyces_cerevisiae': 'TACACAATGTCT'}
STOP_CODONS = ['UGA', 'UAG', 'UAA']
FORBIDDEN_PATTERNS = ["AAAA", "CCCC", "GGGG", "UUUU", "KKKKKK", "MMMMMM", "RRRRRR", "SSSSSS", "WWWWWW", "YYYYYY"]
START_CODON = "AUG"


class SwitchGenerator:
    SEED = 42

    def __init__(self, translated_gene_sequence: str, leader_sequence: str = "GGG",
                 a_domain_size: int = 8, b_domain_size: int = 15, extra_stem_length: int = 0,
                 loop_size: int = 12, linker_pattern: str = "",
                 nucleotides_count_to_trim_gene_after: int = 30):
        self.model = Model(material='rna', celsius=37)
        self.translated_gene_sequence = translated_gene_sequence[:nucleotides_count_to_trim_gene_after]
        self.b_domain_size = b_domain_size
        self.extra_stem_length = extra_stem_length
        self.loop_size = loop_size
        self.leader = Domain(leader_sequence, name="leader", material="rna")
        self.a_domain = Domain(f"N{a_domain_size}", name="a", material="rna")
        self.b_domain = Domain(f"N{b_domain_size}", name="b", material="rna")
        self.extra_stem = Domain(f"N{extra_stem_length}", name="b end", material="rna")
        self.loop = Domain(f"N{loop_size}", name="loop", material="rna")
        self.kozak = Domain(KOZAK_SEQUENCE, name="kozak", material="rna")
        self.linker = Domain(linker_pattern, name="linker", material="rna")
        self.translated_gene = Domain(self.translated_gene_sequence, name="translated gene", material="rna")
        self.trigger_length = a_domain_size + b_domain_size
        self.switch = TargetStrand([
            self.leader, self.a_domain, self.b_domain, self.extra_stem,
            self.loop, ~self.extra_stem, ~self.b_domain, self.kozak, self.linker, self.translated_gene
        ], name="switch")
        self.switch_structure = (f".{len(leader_sequence)}.{a_domain_size}({b_domain_size + extra_stem_length}"
                                 f".{loop_size}){b_domain_size + extra_stem_length}.{len(KOZAK_SEQUENCE)}"
                                 f".{len(linker_pattern)}.{len(self.translated_gene_sequence)}")
        self.switch_complex = TargetComplex([self.switch], self.switch_structure, name='switch complex')
        self.activated_structure = self._get_activated_structure()
        self.stop_codon_pattern = Pattern(STOP_CODONS, scope=[~self.b_domain])
        self.patterns_to_prevent = Pattern(FORBIDDEN_PATTERNS)
        self.options = DesignOptions(f_stop=0.02, seed=self.SEED, wobble_mutations=True)
        self.complex_name = 'activated_complex'

    @property
    def prescribed_stem_structure(self):
        return (f"({self.b_domain_size + self.extra_stem_length}"
                f".{self.loop_size}){self.b_domain_size + self.extra_stem_length}")

    def get_switch(self, trigger_sequence: str, healthy_sequence: str = None) -> Tuple[str, float, float, str]:
        assert len(trigger_sequence) == self.trigger_length

        trigger_domain = Domain(trigger_sequence, name="trigger domain", material='rna')
        trigger = TargetStrand([trigger_domain], name="trigger")
        trigger_structure = "." * len(trigger_sequence)
        trigger_complex = TargetComplex([trigger], trigger_structure, name='trigger complex')

        if healthy_sequence:
            healthy_sequence_domain = Domain(healthy_sequence, name="healthy trigger domain", material='rna')
            healthy_trigger = TargetStrand([healthy_sequence_domain], name="healthy trigger")
            healthy_trigger_structure = "." * len(healthy_sequence)
            healthy_trigger_complex = TargetComplex([healthy_trigger], healthy_trigger_structure,
                                                    name='healthy trigger complex')

        activated_structure = self._get_activated_structure()
        activated_complex = TargetComplex([trigger, self.switch], activated_structure, name=self.complex_name)

        reactants = TargetTube(on_targets={self.switch_complex: 1e-08, trigger_complex: 1e-08},
                               off_targets=SetSpec(max_size=2, exclude=tuple([activated_complex])), name='reactants')

        if healthy_sequence:
            off_target_reactants = TargetTube(on_targets={self.switch_complex: 1e-08, healthy_trigger_complex: 1e-08},
                                              off_targets=SetSpec(max_size=2), name='off_target_reactants')

        products = TargetTube(on_targets={activated_complex: 1e-08}, off_targets=SetSpec(max_size=2), name='products')

        if healthy_sequence:
            tubes = [reactants, products, off_target_reactants]
        else:
            tubes = [reactants, products]

        my_design = tube_design(tubes=tubes,
                                # hard_constraints=[self.stop_codon_pattern, self.patterns_to_prevent],
                                hard_constraints=[self.stop_codon_pattern],
                                soft_constraints=[],
                                model=self.model, options=self.options)

        design_results = my_design.run(trials=1)[0]

        switch_designed_strand = [v for k, v in design_results.to_analysis.strands.items() if
                                  v.name == "switch"].pop()
        ensemble_defect = self._get_ensemble_defect(design_results)
        complex_concentration = self._get_trigger_switch_complex_concentration(design_results, self.complex_name)
        switch_mfe_structure = self._get_switch_mfe_structure(design_results, self.switch_complex)
        return str(switch_designed_strand), ensemble_defect, complex_concentration, switch_mfe_structure

    def _get_activated_structure(self) -> str:
        return ('(' * self.trigger_length + '+' + '.' * self.leader.nt() + ')' * self.trigger_length +
                '.' * (self.switch.nt() - self.leader.nt() - self.trigger_length))

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
