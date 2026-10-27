"""Per-run owners of executable medical-tool resources."""

import os

from mendelprime.domain import require

from . import phenotype as phenotype_tools


class MedicalTools:
    def __init__(self, knowledge, cache_root, external, *, phenotype=False, genotype=None, work_root=None):
        self.knowledge, self.cache_root, self.external = knowledge, cache_root, external
        self._genotype, self.work_root = genotype, work_root
        self._secret = os.environ.get(external["amelie_pseudonym_secret_env"], "")
        self._phenotype = (
            phenotype_tools.load_phenotype_resources(knowledge.resources["resources"]) if phenotype else None
        )

    def profile(self, patient):
        phenotype_tools.require_positive_hpo(patient)
        return phenotype_tools.build_profile(patient, self.knowledge.hpo())

    def phenotype_candidates(self, profile):
        require(self._phenotype is not None, "tools", "phenotype_tools_required")
        return phenotype_tools.collect_candidates(
            self._phenotype, profile, self.cache_root, self.external, self._secret
        )

    def genotype(self, patient):
        from .native import acquire

        require(
            self._genotype is not None and self.work_root is not None,
            "genotype_tools",
            "genotype_tools_required",
        )
        return acquire(patient, self._genotype, self.work_root)
