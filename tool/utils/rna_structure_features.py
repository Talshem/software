

class RNAStructureFeatures(object):
    def __init__(self):
        pass

    @staticmethod
    def dot_parens_plus_fold_ratio(dpp_str):
        n_dots = dpp_str.count('.')
        tot = len(dpp_str)
        fold_ratio = (tot - n_dots) / tot
        return fold_ratio

    