#!/usr/bin/env python
import gc
import numpy as np
import pandas as pd
import pickle
import os
import math
import torch
import warnings
import yaml
import sys
from build_vocab import WordVocab
from utils import split
from transformers import T5EncoderModel, T5Tokenizer
from pretrain_trfm import TrfmSeq2seq

warnings.filterwarnings(action="ignore", category=UserWarning)

inputs_path = ""

if not inputs_path:
    inputs_path = os.getenv("INPUTS")  # From your env.sh file
    if inputs_path is None:
        raise ValueError("The INPUTS environment variable is not set.")

inputs_file = os.path.join(inputs_path, "inputs.yml")
with open(inputs_file, "r") as file:
    data = yaml.safe_load(file)
if not os.path.isfile(inputs_file):
    raise FileNotFoundError(
        f"The 'inputs.yml' file could not be found at {inputs_file}."
    )

# Input and output file paths
output_file_path = os.path.join(inputs_path, data["output_file_path"])
os.makedirs(output_file_path, exist_ok=True)
transporters = data["transporters"]

# UNIKP Python Libraries from bash environment variable
model_path = os.environ.get("UNIKP")
sys.path.append(model_path)

seqs_smiles_df = pd.read_csv(os.path.join(output_file_path, "sequences_smiles.csv"))

with open(os.path.join(model_path, "UniKP for kcat.pkl"), "rb") as f:
    model = pickle.load(f)


# Function to check if a pair of Substrate Smiles and Sequence has been assessed
def is_assessed(substrate_smiles, sequence):
    return (substrate_smiles, sequence) in predicted_pair


def contains_keywords(cell):
    return any(keyword.lower() in str(cell).lower() for keyword in transporters)


# Collect sequences and smiles in batches
sequences = []
smiles = []
indices = []
batch = 0
batch_len = 20
predicted_pair = {}

for index, row in seqs_smiles_df.iterrows():
    if (
        type(row["Sequence"]) is not float
        and not contains_keywords(row["Reaction name"])
        and np.isnan(row["Kcat"])
    ):
        if row["Substrate Smiles"] != "Compound not found":
            pair_key = (row["Substrate Smiles"], row["Sequence"])
            if is_assessed(*pair_key):
                seqs_smiles_df.at[index, "Kcat"] = predicted_pair[pair_key]
            else:
                sequences.append(row["Sequence"])
                smiles.append(row["Substrate Smiles"])
                indices.append(index)


def process_sequence(seq):
    if len(seq) > 1000:
        return seq[:500] + seq[-500:]
    return seq


def Seq_to_vec(Sequence):
    sequences_Example = [" ".join(process_sequence(seq)) for seq in Sequence]
    num_sequences = len(sequences_Example)

    tokenizer = T5Tokenizer.from_pretrained(
        os.path.join(model_path, "prot_t5_xl_uniref50"), do_lower_case=False
    )
    model = T5EncoderModel.from_pretrained(
        os.path.join(model_path, "prot_t5_xl_uniref50")
    )
    gc.collect()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if torch.cuda.is_available():
        print("Let's use", torch.cuda.device_count(), "GPUs!")
        model = torch.nn.DataParallel(model)
        batch_size = min(
            num_sequences, torch.cuda.device_count() * 8
        )  # Adjust batch size for multiple GPUs
    else:
        print("Let's use", os.cpu_count(), "CPUs!")
        batch_size = min(num_sequences, os.cpu_count())
    model = model.to(device).eval()

    features = []

    for i in range(0, num_sequences, batch_size):
        batch_sequences = sequences_Example[i : i + batch_size]
        batch_ids = tokenizer.batch_encode_plus(
            batch_sequences, add_special_tokens=True, padding=True
        )
        input_ids = torch.tensor(batch_ids["input_ids"]).to(device)
        attention_mask = torch.tensor(batch_ids["attention_mask"]).to(device)

        with torch.no_grad():
            embedding = model(input_ids=input_ids, attention_mask=attention_mask)

        embedding = embedding.last_hidden_state
        for seq_num in range(embedding.size(0)):
            seq_len = (attention_mask[seq_num] == 1).sum()
            seq_emd = embedding[seq_num][: seq_len - 1]
            features.append(seq_emd)

    print("Finished for sequence tokenizer loop")

    return features


# Process all sequences and smiles in one go
if sequences:
    features = Seq_to_vec(sequences)


def normalize_feature(features):
    # Perform normalization directly on the GPU
    features_normalize = torch.stack([f.mean(dim=0) for f in features], dim=0)

    # Move features_normalize back to CPU if needed
    features_normalize = features_normalize.cpu().numpy()
    return features_normalize


if features:
    seq_vecs = normalize_feature(features)


def smiles_to_vec(Smiles):
    pad_index = 0
    unk_index = 1
    eos_index = 2
    sos_index = 3
    mask_index = 4
    vocab = WordVocab.load_vocab(os.path.join(model_path, "vocab.pkl"))
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    def get_inputs(sm):
        seq_len = 220
        sm = sm.split()
        if len(sm) > 218:
            print("SMILES is too long ({:d})".format(len(sm)))
            sm = sm[:109] + sm[-109:]
        ids = [vocab.stoi.get(token, unk_index) for token in sm]
        ids = [sos_index] + ids + [eos_index]
        seg = [1] * len(ids)
        padding = [pad_index] * (seq_len - len(ids))
        ids.extend(padding), seg.extend(padding)
        return ids, seg

    def get_array(smiles):
        x_id, x_seg = [], []
        for sm in smiles:
            a, b = get_inputs(sm)
            x_id.append(a)
            x_seg.append(b)
        return torch.tensor(x_id).to(device), torch.tensor(x_seg).to(device)

    trfm = TrfmSeq2seq(len(vocab), 256, len(vocab), 4)
    # Use map_location to ensure the model is loaded on the device if GPU is not available
    trfm.load_state_dict(
        torch.load(os.path.join(model_path, "trfm_12_23000.pkl"), map_location=device)
    )
    trfm.to(device)
    trfm.eval()

    X = []
    for smile in Smiles:
        x_split = [split(smile)]
        xid, xseg = get_array(x_split)
        X.append(trfm.encode(torch.t(xid).to(device))[0])
    return X


if smiles:
    smiles_vecs = smiles_to_vec(smiles)

if sequences and smiles:
    fused_vectors = np.concatenate((smiles_vecs, seq_vecs), axis=1)
    pre_kcats = model.predict(fused_vectors)
    kcates = [math.pow(10, pre_kcats[i]) for i in range(len(pre_kcats))]

    for i, index in enumerate(indices):
        seqs_smiles_df.at[index, "Kcat"] = kcates[i]
        pair_key = (smiles[i], sequences[i])
        predicted_pair[pair_key] = kcates[i]

    # Save the DataFrame periodically
    batch += 1
    if batch == batch_len:
        seqs_smiles_df.to_csv(
            os.path.join(output_file_path, "sequences_smiles_complete.csv"), index=False
        )
        batch = 0

seqs_smiles_df.to_csv(
    os.path.join(output_file_path, "sequences_smiles_complete.csv"), index=False
)
