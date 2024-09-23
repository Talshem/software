
# # UFOLD
# # https://ufold.ics.uci.edu/
# import math
# import numpy as np
# import time
# from datetime import timedelta
# from utils.enegry_calculator import EnergyCalculator
# from utils.rna_structure_features import RNAStructureFeatures
# import utils.seqeunce_consts as sc

# # Didi code for trigger selection from Trigate project


# class RNATriggerSelector(object):
#     TRIGGER_SCORE_COL_NAME = 'trigger_score'

#     DEFAULT_SLIDE_WINDOW_SIZE = 50
#     DEFAULT_TRIGGER_SIZE = 36
#     DEFAULT_NEIGHBOURHOOD_WIN_SIZE = 11

#     def __init__(
#             self,
#             seq,
#             aux_seq=None,
#             win_size=DEFAULT_SLIDE_WINDOW_SIZE,
#             trigger_size=DEFAULT_TRIGGER_SIZE,
#             neighbourhood_win_size=DEFAULT_NEIGHBOURHOOD_WIN_SIZE,):

#         self.seq = seq
#         self.aux_seq = aux_seq

#         self.win_size = win_size
#         self.trigger_size = trigger_size
#         self.neighbourhood_win_size = neighbourhood_win_size

#         self.triggers_avg_fold_energies = None
#         self.triggers_avg_avg_fold_energies = None
#         self.triggers_avg_fold_prob = None

#         self.aux_energies = None
#         self.aux_probs = None

#         # for debug view
#         self.core_scores = None
#         self.aux_scores = None

#     def calculate_single_rna_triggers_score(self):
#         self.calculate_single_rna_triggers_score_core()


#     def calculate_single_rna_triggers_score_core(self):

#         print('calculate_single_rna_triggers_score_core')
#         start_time = time.monotonic()
#         seq = self.seq
#         win_size = self.win_size
#         trigger_size = self.trigger_size
#         seq_size = len(seq)


#         neighbourhood_win_size = self.neighbourhood_win_size

#         assert win_size >= trigger_size
#         assert win_size >= neighbourhood_win_size
#         assert neighbourhood_win_size % 2 == 1

#         n_triggers = seq_size - trigger_size + 1
#         n_measurements = win_size - trigger_size + 1

#         triggers_fold_energies = np.zeros((n_triggers, n_measurements), dtype=np.float32)
#         triggers_mfe_structure_nuc_fold_ratio = np.zeros((n_triggers, n_measurements), dtype=np.float32)
#         triggers_ensemble_nuc_fold_avg_prob = np.zeros((n_triggers, n_measurements), dtype=np.float32)
#         triggers_n_measurements = np.zeros(n_triggers, dtype=np.int32)

#         assert len(seq) > win_size

#         n_windows = len(seq) - win_size + 1


#         for wind_ind in range(n_windows):
#             part_seq_start = wind_ind
#             part_seq_end = wind_ind + win_size
#             assert part_seq_start >= 0
#             assert part_seq_end <= len(seq)


#             part_seq = seq[part_seq_start: part_seq_end]

#             results = EnergyCalculator.calc_self_mfe_structure_and_fold_energy_and_probabilities(part_seq)

#             free_energy = results.free_energy
#             # mfe = results.mfe[0].energy

#             for inner_win_ind in range(n_measurements):
#                 trigger_ind = wind_ind + inner_win_ind
#                 current_repeat = triggers_n_measurements[trigger_ind]

#                 triggers_fold_energies[trigger_ind, current_repeat] = free_energy

#                 inner_trigger_start_pos = inner_win_ind
#                 inner_trigger_end_pos = inner_win_ind + trigger_size
#                 trigger_sub_struct = str(results.mfe[0].structure)[inner_trigger_start_pos:inner_trigger_end_pos]

#                 triggers_mfe_structure_nuc_fold_ratio[trigger_ind, current_repeat] = \
#                     RNAStructureFeatures.dot_parens_plus_fold_ratio(trigger_sub_struct)

#                 avg_fold_prob = np.mean(results.pairs.array.diagonal()[inner_trigger_start_pos:inner_trigger_end_pos])
#                 triggers_ensemble_nuc_fold_avg_prob[trigger_ind, current_repeat] = avg_fold_prob

#                 # TODO CG ratio of folded VS CG ratio folded in trigger
#                 # TODO calculate the only trigger sub strand with other parts it folds with as sub strands
#                 #  and calculate co-fold energy more accurate than taking partial of window energy

#                 triggers_n_measurements[trigger_ind] += 1

#         triggers_fold_energies = \
#             triggers_fold_energies * \
#             triggers_mfe_structure_nuc_fold_ratio * \
#             triggers_ensemble_nuc_fold_avg_prob

#         triggers_avg_fold_prob = np.zeros(n_triggers, dtype=np.float32)
#         triggers_avg_fold_energies = np.zeros(n_triggers, dtype=np.float32)
#         for trigger_ind in range(n_triggers):
#             n_measurements = triggers_n_measurements[trigger_ind]

#             vals = triggers_fold_energies[trigger_ind, :n_measurements]
#             triggers_avg_fold_energies[trigger_ind] = np.mean(vals)

#             vals = triggers_ensemble_nuc_fold_avg_prob[trigger_ind, :n_measurements]
#             triggers_avg_fold_prob[trigger_ind] = np.mean(vals)

#         pad_size = (neighbourhood_win_size - 1) // 2
#         start_pad_val = triggers_avg_fold_energies[0]
#         end_pad_val = triggers_avg_fold_energies[-1]
#         padded_triggers_avg_fold_energies = \
#             ([start_pad_val] * pad_size) + \
#             triggers_avg_fold_energies.tolist() + \
#             ([end_pad_val] * pad_size)


#         # import scipy.ndimage.filters as ndif
#         # triggers_neighbourhood_avg_fold_energies = ndif.uniform_filter1d(
#         #     padded_triggers_avg_fold_energies, neighbourhood_win_size,
#         #     mode='constant',
#         #     origin=-(neighbourhood_win_size // 2))[:-(neighbourhood_win_size - 1)]

#         triggers_neighbourhood_avg_fold_energies = np.convolve(
#             padded_triggers_avg_fold_energies,
#             np.ones(neighbourhood_win_size) / float(neighbourhood_win_size), 'valid')

#         self.triggers_avg_fold_energies = triggers_avg_fold_energies
#         self.triggers_avg_avg_fold_energies = triggers_neighbourhood_avg_fold_energies
#         self.triggers_avg_fold_prob = triggers_avg_fold_prob

#         end_time = time.monotonic()
#         print(timedelta(seconds=end_time - start_time))


#     def generate_trigger_score_table(self, is_energy_score=True, n_best=10, file_name=None):
#         seq = self.seq
#         core_scores = 1.0 - self.triggers_avg_fold_prob
#         if is_energy_score:
#             core_scores = 1.0 / np.abs(1e-3 - self.triggers_avg_fold_energies)
#             # core_scores = self.triggers_avg_avg_fold_energies
#             # core_scores = NormUtil.min_max_scaler(core_scores)

#         self.core_scores = core_scores

#         total_scores = core_scores

#         aux_scores = None
#         if self.aux_seq:
#             aux_scores = 1.0 - self.aux_probs
#             if is_energy_score:
#                 aux_scores = 1.0 / np.abs(1e-3 - self.aux_energies)
#             # norm_aux_score = NormUtil.min_max_scaler(aux_score)

#             # for usage in view
#             self.aux_scores = aux_scores

#             total_scores = total_scores * aux_scores


#         # TODO understand why NUMPY BUG += *= !!!!
#         # dummy = self.norm_fold_energies
#         # dummy += self.norm_aux_score
#         # a = dummy[0]
#         # b = self.norm_fold_energies[0]
#         # c = self.norm_fold_energies[0]
#         # d = self.norm_aux_score[0]
#         # assert a != b

#         # indexes = np.argsort(total_scores)[::-1][:n_best]
#         indexes = np.argsort(total_scores)[::-1]

#         dist_threshold = max(5, len(self.seq) // (2 * n_best), self.trigger_size // 5)

#         arr = np.array([], dtype=int)
#         arr = np.append(arr, indexes[0])

#         best_indexes = [indexes[0]]

#         for ind in indexes:
#             near_ind = self.find_nearest(arr, ind)
#             abs_dist = math.fabs(ind - near_ind)
#             if abs_dist < dist_threshold:
#                 continue
#             best_indexes.append(ind)
#             arr = np.append(arr, ind)
#             arr.sort()
#             if arr.shape == n_best:
#                 break

#         indexes = best_indexes
#         best_triggers = [seq[i:i+self.trigger_size] for i in indexes]
#         best_aux_scores = None
#         if self.aux_seq:
#             best_aux_scores = aux_scores[indexes]

#         best_core_scores = core_scores[indexes]
#         best_total_scores = total_scores[indexes]

#         d = {
#             'start_index:': indexes,
#             'trigger': best_triggers,
#             self.TRIGGER_SCORE_COL_NAME: best_total_scores,
#             'trigger_core_scores': best_core_scores
#         }

#         if self.aux_seq:
#             d.update({'trigger_aux_score': best_aux_scores})

#         import pandas as pd

#         df = pd.DataFrame.from_dict(d)

#         if file_name:
#             df.to_csv(path_or_buf=file_name, index=True)

#         return df

#         # import os
#         # dir_path = os.path.dirname(os.path.realpath(__file__))
#         # i = 3

#     def view_results(self, file_path=None):
#         # energies = self.triggers_avg_fold_energies
#         # avg_energies = self.triggers_avg_avg_fold_energies
#         # probs = self.triggers_avg_fold_prob
#         # aux = self.aux_score

#         # assert len(energies) == len(probs)

#         import matplotlib.pyplot as plt

#         plt.rcParams['figure.dpi'] = 600

#         n_plots = 1
#         if self.aux_seq:
#             n_plots += 1

#         fig, axes = plt.subplots(1, n_plots)

#         fig.suptitle('trigger scores')

#         x = range(len(self.core_scores))

#         ax_no = 0
#         ax_no = self.add_plot(axes, n_plots, ax_no, x, self.core_scores, 'core scores', 'offset', 'score')
#         # axes[ax_no].plot(x, energies, label='energies')
#         # axes[ax_no].set_title('energies')
#         # ax_no += 1

#         # axes[ax_no].plot(x, avg_energies, label='avg energies')
#         # axes[ax_no].set_title('avg energies')

#         # ax_no = self.add_plot(axes, ax_no, x, probs, 'probs')
#         # axes[2].plot(x, probs, label='probs')
#         # axes[2].set_title('probs')

#         if self.aux_seq:
#             ax_no = self.add_plot(axes, n_plots, ax_no, x, self.aux_scores, 'aux scores', 'offset', 'score')
#             # axes[3].plot(x, aux, label='aux')
#             # axes[3].set_title('aux')

#         # plt.xlabel('trigger pos')
#         # plt.ylabel('value')
#         # plt.title('graph')

#         plt.legend()

#         if file_path:
#             plt.savefig("dump.jpg", dpi=600)
#         else:
#             plt.show()

#     @staticmethod
#     def add_plot(axes, n_plots, ax_no, x, y, title, x_label, y_label):
#         if n_plots > 1:
#             axes[ax_no].plot(x, y, label=title)
#             axes[ax_no].set_title(title)
#             axes[ax_no].set(xlabel=x_label, ylabel=y_label)
#         else:
#             axes.plot(x, y, label=title)
#             axes.set_title(title)
#             axes.set(xlabel=x_label, ylabel=y_label)
#         return ax_no + 1

#     @staticmethod
#     def find_nearest(array, value):
#         idx = np.searchsorted(array, value, side="left")
#         if idx > 0 and (idx == len(array) or math.fabs(value - array[idx - 1]) < math.fabs(value - array[idx])):
#             return array[idx - 1]
#         else:
#             return array[idx]


# if __name__ == '__main__':

#     win_len_ = 120
#     nei_len_ = 11
#     trigger_len_ = 30

#     out_file_name = f"triggers_w{win_len_}_n{nei_len_}_t{trigger_len_}.csv"


#     m_cherry_orf = sc.M_CHERRY_ORF
#     rts_ = RNATriggerSelector(win_size=win_len_,
#                               trigger_size=trigger_len_,
#                               neighbourhood_win_size=nei_len_,
#                               seq=m_cherry_orf,
#                             )
#     rts_.calculate_single_rna_triggers_score()
#     df_ = rts_.generate_trigger_score_table(is_energy_score=True, n_best=10, file_name=out_file_name)

