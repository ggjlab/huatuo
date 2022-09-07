import argparse
import math
import pyfasta
import torch
from torch import nn
import numpy as np
import pandas as pd
import h5py

parser = argparse.ArgumentParser(description='Predict variant chromatin effects')
parser.add_argument('inputfile', type=str, help='Input file in vcf format')
parser.add_argument('--maxshift', action="store",
                    dest="maxshift", type=int, default=800,
                    help='Maximum shift distance for computing nearby effects')
parser.add_argument('--inputsize', action="store", dest="inputsize", type=int,
                    default=2000, help="The input sequence window size for neural network")
parser.add_argument('--batchsize', action="store", dest="batchsize",
                    type=int, default=32, help="Batch size for neural network predictions.")
parser.add_argument('--cuda', action='store_true')
parser.add_argument('--hg19_fa', action="store",
                    dest="hg19_fa", default='./resources/hg19.fa')
parser.add_argument('--deepsea', action="store",
                    dest="deepsea", default='./resources/deepsea.beluga.pth')

args = parser.parse_args()

genome = pyfasta.Fasta(args.hg19_fa)

class LambdaBase(nn.Sequential):
    def __init__(self, fn, *args):
        super(LambdaBase, self).__init__(*args)
        self.lambda_func = fn

    def forward_prepare(self, input):
        output = []
        for module in self._modules.values():
            output.append(module(input))
        return output if output else input

class Lambda(LambdaBase):
    def forward(self, input):
        return self.lambda_func(self.forward_prepare(input))

class Beluga(nn.Module):
    def __init__(self):
        super(Beluga, self).__init__()
        self.model = nn.Sequential(
            nn.Sequential(
                nn.Conv2d(4,320,(1, 8)),
                nn.ReLU(),
                nn.Conv2d(320,320,(1, 8)),
                nn.ReLU(),
                nn.Dropout(0.2),
                nn.MaxPool2d((1, 4),(1, 4)),
                nn.Conv2d(320,480,(1, 8)),
                nn.ReLU(),
                nn.Conv2d(480,480,(1, 8)),
                nn.ReLU(),
                nn.Dropout(0.2),
                nn.MaxPool2d((1, 4),(1, 4)),
                nn.Conv2d(480,640,(1, 8)),
                nn.ReLU(),
                nn.Conv2d(640,640,(1, 8)),
                nn.ReLU(),
            ),
            nn.Sequential(
                nn.Dropout(0.5),
                Lambda(lambda x: x.view(x.size(0),-1)),
                nn.Sequential(Lambda(lambda x: x.view(1,-1) if 1==len(x.size()) else x ),nn.Linear(67840,2003)),
                nn.ReLU(),
                nn.Sequential(Lambda(lambda x: x.view(1,-1) if 1==len(x.size()) else x ),nn.Linear(2003,2002)),
            ),
            nn.Sigmoid(),
        )

    def forward(self, x):
        return self.model(x)

model = Beluga()
model.load_state_dict(torch.load(args.deepsea))
model.eval()
if args.cuda:
    model.cuda()

CHRS = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
        'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
        'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX','chrY']


inputfile = args.inputfile
maxshift = args.maxshift
inputsize = args.inputsize
batchSize = args.batchsize
windowsize = inputsize + 100


def encodeSeqs(seqs, inputsize=2000):

    seqsnp = np.zeros((len(seqs), 4, inputsize), np.bool_)

    mydict = {'A': np.asarray([1, 0, 0, 0]), 'G': np.asarray([0, 1, 0, 0]),
            'C': np.asarray([0, 0, 1, 0]), 'T': np.asarray([0, 0, 0, 1]),
            'N': np.asarray([0, 0, 0, 0]), 'H': np.asarray([0, 0, 0, 0]),
            'a': np.asarray([1, 0, 0, 0]), 'g': np.asarray([0, 1, 0, 0]),
            'c': np.asarray([0, 0, 1, 0]), 't': np.asarray([0, 0, 0, 1]),
            'n': np.asarray([0, 0, 0, 0]), '-': np.asarray([0, 0, 0, 0])}

    n = 0
    for line in seqs:
        cline = line[int(math.floor(((len(line) - inputsize) / 2.0))):int(math.floor(len(line) - (len(line) - inputsize) / 2.0))]
        for i, c in enumerate(cline):
            seqsnp[n, :, i] = mydict[c]
        n = n + 1

    # get the complementary sequences
    dataflip = seqsnp[:, ::-1, ::-1]
    seqsnp = np.concatenate([seqsnp, dataflip], axis=0)
    return seqsnp


def fetchSeqs(chr, pos, ref, alt, shift=0, inputsize=2000):

    windowsize = inputsize + 100
    mutpos = int(windowsize / 2 - 1 - shift)
    # return string: ref sequence, string: alt sequence, Bool: whether ref allele matches with reference genome
    seq = genome.sequence({'chr': chr, 'start': pos + shift -
                           int(windowsize / 2 - 1), 'stop': pos + shift + int(windowsize / 2)})
    return seq[:mutpos] + ref + seq[(mutpos + len(ref)):], seq[:mutpos] + alt + seq[(mutpos + len(ref)):], seq[mutpos:(mutpos + len(ref))].upper() == ref.upper()


vcf = pd.read_csv(inputfile, sep='\t', header=None, comment='#')

# standardize
vcf.iloc[:, 0] = 'chr' + vcf.iloc[:, 0].map(str).str.replace('chr', '')
vcf = vcf[vcf.iloc[:, 0].isin(CHRS)]

for shift in [0, ] + list(range(-200, -maxshift - 1, -200)) + list(range(200, maxshift + 1, 200)):
    refseqs = []
    altseqs = []
    ref_matched_bools = []
    for i in range(vcf.shape[0]):
        refseq, altseq, ref_matched_bool = fetchSeqs(
            vcf[0][i], vcf[1][i], vcf[2][i], vcf[3][i], shift=shift, inputsize=inputsize)
        refseqs.append(refseq)
        altseqs.append(altseq)
        ref_matched_bools.append(ref_matched_bool)

    if shift == 0:
        # only need to be checked once
        print("Number of variants with reference allele matched with reference genome:")
        print(np.sum(ref_matched_bools))
        print("Number of input variants:")
        print(len(ref_matched_bools))

    ref_encoded = encodeSeqs(refseqs, inputsize=inputsize).astype(np.float32)
    alt_encoded = encodeSeqs(altseqs, inputsize=inputsize).astype(np.float32)

    ref_preds = []
    for i in range(int(1 + (ref_encoded.shape[0]-1) / batchSize)):
        input = torch.from_numpy(ref_encoded[int(i*batchSize):int((i+1)*batchSize),:,:]).unsqueeze(2)
        if args.cuda:
            input = input.cuda()
        ref_preds.append(model.forward(input).cpu().detach().numpy().copy())
    ref_preds = np.vstack(ref_preds)

    alt_preds = []
    for i in range(int(1 + (alt_encoded.shape[0]-1) / batchSize)):
        input = torch.from_numpy(alt_encoded[int(i*batchSize):int((i+1)*batchSize),:,:]).unsqueeze(2)
        if args.cuda:
            input = input.cuda()
        alt_preds.append(model.forward(input).cpu().detach().numpy().copy())
    alt_preds = np.vstack(alt_preds)

    diff = alt_preds - ref_preds
    f = h5py.File(inputfile + '.shift_' + str(shift) + '.diff.h5', 'w')
    f.create_dataset('pred', data=diff)
    f.close()
