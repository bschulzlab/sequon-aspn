import pandas as pd
import re
from get_uniprot import UniprotSequence, UniprotParser
from io import StringIO
from sequence import Sequence, sequon_re, sequon_re_modded, Peptide

accession_re = re.compile(
    "(?P<accession>[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(?P<isotype>-\d)?")

normalized_result = pd.read_csv("data/Normalised_result_AspN_No_Reverse.txt", sep="\t")

acc_list = []

for i in normalized_result["Protein"].unique():
    acc = UniprotSequence(i, True)
    if acc.accession != "":
        acc_list.append(UniprotSequence(i, True))

p = UniprotParser(acc_list)
data = []
for i in p.parse():
    seq = {}
    for line in StringIO(i):
        if line[0] == ">":
            if seq:
                data.append(seq)
            seq = {}
            id_full = line[1: line.find(" ")]
            id_split = id_full.split("|")
            seq["id"] = id_split[1]
            seq["entry name"] = id_split[2]
            seq["full id"] = id_full
            seq["sequence"] = ""
        else:
            if seq:
                seq["sequence"] += line.strip()
    data.append(seq)

data = pd.DataFrame(data)
data = data.set_index("full id")
# normalized_result = pd.read_csv("data/Normalised_result_AspN_No_Reverse.txt", sep="\t")

for i, df in normalized_result.groupby("Protein"):
    if i in data.index:
        print(i)
        whole_sequon_data = []
        seq = Sequence(data.loc[i]["sequence"])
        seq_length = len(seq.seq)
        for match in seq.find_with_regex(sequon_re, seq.gaps()):
            whole_sequon_data.append([match.start, match.stop, seq.seq[match.start: match.stop], "", "", ""])
        whole_sequon_data = pd.DataFrame(whole_sequon_data, columns=[
            "Start", "Stop", "Sequence", "Start_modded", "Stop_modded", "Modification"])
        #print(whole_sequon_data)
        for ind, r in df.iterrows():
            s = Peptide(r["Peptide"])

            original_position, stop_position, extra_seq = s.map_seq(seq.seq)
            sub_seq = Sequence(s.seq + extra_seq)
            sequon = whole_sequon_data[(whole_sequon_data["Start"] >= original_position) &
                                       (whole_sequon_data["Start"] < stop_position)]

            if not sequon.empty:
                modded_sequon = []
                #print(sequon)
                current_sequon = -1
                gap = 0

                for match in sub_seq.find_with_regex(sequon_re_modded, sub_seq.gaps()):

                    if sub_seq.seq[match.start] == "[":

                        if current_sequon != -1:
                            remapped = original_position + match.start - gap - 1
                            #print(remapped)
                            if remapped == sequon.iloc[current_sequon]["Start"]:
                                sequon.iloc[current_sequon, 5] = sub_seq.seq[match.start + 1: match.stop - 1]
                            gap += match.stop - match.start
                    else:
                        current_sequon += 1
            else:
                # print(stop_position)
                if stop_position in whole_sequon_data["Start"].values:
                    # print(whole_sequon_data[whole_sequon_data["Start"] == stop_position])
                    a = whole_sequon_data[whole_sequon_data["Start"] == stop_position]
                    normalized_result.at[ind, "Sequon After C-term"] = a["Sequence"].values[0] + "(" + str(stop_position+1) + ")"
            #print(sequon)
            for n, n_sequon in sequon.iterrows():
                if "Sequon Position" not in normalized_result.columns:
                    normalized_result.at[ind, "Sequon Position"] = str(
                        n_sequon["Start"]+1) + ";"

                    normalized_result.at[ind, "Sequon Sequence"] = n_sequon["Sequence"] \
                                                            + "(" + str(n_sequon["Start"]+1) + ")" + ";"
                elif pd.isnull(normalized_result.at[ind, "Sequon Position"]):
                    normalized_result.at[ind, "Sequon Position"] = str(
                        n_sequon["Start"]+1) + ";"

                    normalized_result.at[ind, "Sequon Sequence"] = n_sequon["Sequence"] \
                                                            + "(" + str(n_sequon["Start"]+1) + ")" + ";"
                else:
                    normalized_result.at[ind, "Sequon Position"] = normalized_result.at[ind, "Sequon Position"] + str(
                        n_sequon["Start"]+1) + ";"
                    normalized_result.at[ind, "Sequon Sequence"] = normalized_result.at[ind, "Sequon Sequence"] + n_sequon["Sequence"] \
                                                            + "(" + str(n_sequon["Start"]+1) + ")" + ";"
                if "Sequon Modification" not in normalized_result.columns:
                    normalized_result.at[ind, "Sequon Modification"] = ""
                if n_sequon["Modification"] != "" and pd.notnull(n_sequon["Modification"]):
                    if "Sequon Modification" not in normalized_result.columns:
                        normalized_result.at[ind, "Sequon Modification"] = n_sequon[
                                                                               "Modification"] + "(" + str(
                            n_sequon["Start"]+1) + ")" + ";"
                    elif pd.isnull(normalized_result.at[ind, "Sequon Modification"]):
                        normalized_result.at[ind, "Sequon Modification"] = n_sequon[
                                                                               "Modification"] + "(" + str(
                            n_sequon["Start"]+1) + ")" + ";"
                    else:
                        normalized_result.at[ind, "Sequon Modification"] = normalized_result.at[
                                                                               ind, "Sequon Modification"] + n_sequon[
                                                                               "Modification"] + "(" + str(
                            n_sequon["Start"]+1) + ")" + ";"

normalized_result.to_csv("result_aspn.txt", sep="\t", index=False)
