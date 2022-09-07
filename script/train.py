import h5py
import argparse
import xgboost as xgb
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler


parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--targetIndex', action="store",
                    dest="targetIndex", type=int)
parser.add_argument('--output', action="store", dest="output")
parser.add_argument('--expFile', action="store", dest="expFile")
parser.add_argument('--inputFile', action="store",
                    dest="inputFile", default='./features.npz')
parser.add_argument('--annoFile', action="store",
                    dest="annoFile", default='./geneanno.csv')
parser.add_argument('--evalFile', action="store",
                     dest="evalFile", default='',help='specify to save holdout set predictions')
parser.add_argument('--filterStr', action="store",
                    dest="filterStr", type=str, default="all")
parser.add_argument('--pseudocount', action="store",
                    dest="pseudocount", type=float, default=0.0001)
parser.add_argument('--num_round', action="store",
                    dest="num_round", type=int, default=40)
parser.add_argument('--l2', action="store", dest="l2", type=float, default=100)
parser.add_argument('--l1', action="store", dest="l1", type=float, default=0)
parser.add_argument('--eta', action="store", dest="eta",
                    type=float, default=.5)
parser.add_argument('--base_score', action="store",
                    dest="base_score", type=float, default=2)
parser.add_argument('--threads', action="store",
                    dest="threads", type=int, default=16)

args = parser.parse_args()

# read resources

geneused = pd.read_csv("./geneused.csv", index_col=0)
geneexp = pd.read_csv(args.expFile,index_col=0)
geneexp = geneexp[~geneexp.index.duplicated(keep='first')]

res = geneexp.index.intersection(geneused.index)

geneexp = geneexp.loc[res]
index = [geneused.index.tolist().index(i) for i in res.tolist()] #
geneused = geneused.loc[res]

geneanno = pd.read_csv(args.annoFile)
geneanno = geneused.merge(geneanno, left_index=True, right_on="symbol", how='left').drop_duplicates(['symbol'])

Xreducedall = np.load(args.inputFile)['features'] 
Xreducedall = Xreducedall[index]
Xreducedall = StandardScaler().fit_transform(Xreducedall)

##
if args.filterStr == 'pc':
    filt = np.asarray(geneanno.iloc[:, -1] == 'protein_coding')
elif args.filterStr == 'lincRNA':
    filt = np.asarray(geneanno.iloc[:, -1] == 'lincRNA')
elif args.filterStr == 'all':
    filt = np.asarray(geneanno.iloc[:, -1] != 'rRNA')
else:
    raise ValueError('filterStr has to be one of all, pc, and lincRNA')

filt = filt * \
    np.isfinite(np.asarray(
        np.log(geneexp.iloc[:, args.targetIndex] + args.pseudocount)))


filt_1 = np.asarray(geneanno.iloc[:, -1] == 'protein_coding')
filt_1 = filt_1 * \
    np.isfinite(np.asarray(
        np.log(geneexp.iloc[:, args.targetIndex] + args.pseudocount)))


# training

trainind = np.asarray(geneanno['seqnames'] != 'chrX') * np.asarray(
    geneanno['seqnames'] != 'chrY') * np.asarray(geneanno['seqnames'] != 'chr8')
testind = np.asarray(geneanno['seqnames'] == 'chr8')

dtrain = xgb.DMatrix(Xreducedall[trainind * filt_1, :])
dtest = xgb.DMatrix(Xreducedall[(testind) * filt, :])

dtrain.set_label(np.asarray(
    np.log(geneexp.iloc[trainind * filt_1, args.targetIndex] + args.pseudocount)))
dtest.set_label(np.asarray(
    np.log(geneexp.iloc[(testind) * filt, args.targetIndex] + args.pseudocount)))

param = {'booster': 'gblinear', 'base_score': args.base_score, 'alpha': 0,
         'lambda': args.l2, 'eta': args.eta, 'objective': 'reg:linear',
         'nthread': args.threads, "early_stopping_rounds": 10}

evallist = [(dtest, 'eval'), (dtrain, 'train')]
num_round = args.num_round
bst = xgb.train(param, dtrain, num_round, evallist)
ypred = bst.predict(dtest)
print(spearmanr(ypred, np.asarray(
     np.log(geneexp.iloc[(testind) * filt, args.targetIndex] + args.pseudocount))))
if args.evalFile != '':
    evaldf = pd.DataFrame({'pred':ypred,'target':np.asarray(
     np.log(geneexp.iloc[(testind) * filt, args.targetIndex] + args.pseudocount))})
    evaldf.to_csv(args.evalFile)
bst.save_model(args.output + args.filterStr + '.pseudocount' + str(args.pseudocount) + '.lambda' + str(args.l2) + '.round' +
               str(args.num_round) + '.basescore' + str(args.base_score) + '.' + geneexp.columns[args.targetIndex] + '.save')
bst.dump_model(args.output + args.filterStr + '.pseudocount' + str(args.pseudocount) + '.lambda' + str(args.l2) + '.round' +
               str(args.num_round) + '.basescore' + str(args.base_score) + '.' + geneexp.columns[args.targetIndex] + '.dump')
