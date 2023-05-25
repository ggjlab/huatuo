import argparse
import xgboost as xgb
import pandas as pd
import numpy as np
import h5py
from six.moves import reduce

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--coorFile', action="store", dest="coorFile")
parser.add_argument('--snpEffectFilePattern', action="store", dest="snpEffectFilePattern",
                    help="SNP effect hdf5 filename pattern. Use SHIFT as placeholder for shifts.")
parser.add_argument('--modelList', action="store", dest="modelList",
                    help="A list of paths of binary xgboost model files (if end with .list) or a combined model file (if ends with .csv).")
parser.add_argument('--nfeatures', action="store",
                    dest="nfeatures", type=int, default=2002)
parser.add_argument('--output', action="store", dest="output")
parser.add_argument('--fixeddist', action="store",
                    dest="fixeddist", default=0, type=int)
parser.add_argument('--maxshift', action="store",
                    dest="maxshift", type=int, default=800)
parser.add_argument('--batchSize', action="store",
                    dest="batchSize", type=int, default=500)
parser.add_argument('--splitFlag', action="store_true", default=False)
parser.add_argument('--splitIndex', action="store",
                    dest="splitIndex", type=int, default=0)
parser.add_argument('--splitFold', action="store",
                    dest="splitFold", type=int, default=10)
parser.add_argument('--threads', action="store", dest="threads",
                    type=int, default=16, help="Number of threads.")
args = parser.parse_args()


def compute_effects(snpeffects, n_snps, all_models, maxshift=800, nfeatures=2002, batchSize=500, old_format=False):

    Xreducedall_diffs = []
    dists = [0, ] + list(range(-200, -maxshift - 1, -200)) + list(range(200, maxshift + 1, 200))

    effect = np.zeros((n_snps, len(all_models))) # (n_snps, n_models)

    for i in range(int( (n_snps - 1) / batchSize) + 1):
        print("Processing " + str(i) + "th batch of "+str(batchSize))
        # compute gene expression change with models
        diff = reduce(lambda x, y: x + y,
                        [np.tile(np.asarray(snpeffects[j][i * batchSize:(i + 1) * batchSize, :]), 191) for j in range(len(dists))] # (n_dists, n_snps, 10)
                        ) # (n_snps, n_features*10)
        if old_format:
            # backward compatibility
            diff = np.concatenate([np.zeros((diff.shape[0], 10, 1)), diff.reshape(
                (-1, 10, 2002))], axis=2).reshape((-1, 20030))
        dtest_ref = xgb.DMatrix(diff * 0)
        dtest_alt = xgb.DMatrix(diff)

        for j in range(len(all_models)):
            effect[i * batchSize:(i + 1) * batchSize, j] = all_models[j].predict(dtest_alt) - \
                            all_models[j].predict(dtest_ref)

    return effect


#load resources
modelList = pd.read_csv(args.modelList,sep='\t',header=0)
models = []
for file in modelList['ModelName']:
        bst = xgb.Booster({'nthread': args.threads})
        bst.load_model(file.strip())
        models.append(bst)

# backward compatibility with earlier model format
if len(bst.get_dump()[0].split('\n')) == 20034:
    old_format = True
else:
    old_format = False


#load input data
maxshift = int(args.maxshift)
snpEffects = []
for shift in [str(n) for n in [0, ] + list(range(-200, -maxshift - 1, -200)) + list(range(200, maxshift + 1, 200))]:
    h5f = h5py.File(args.snpEffectFilePattern.replace(
        'SHIFT', shift), 'r')['/pred']

    if args.splitFlag:
        index_start = int((args.splitIndex - 1) *
                          np.ceil(float(h5f.shape[0] / 2) / args.splitFold))
        index_end = int(np.minimum(
            (args.splitIndex) * np.ceil(float(h5f.shape[0] / 2) / args.splitFold), (h5f.shape[0] / 2)))
    else:
        index_start = 0
        index_end = int(h5f.shape[0] / 2)

    snp_temp = (np.asarray(h5f[index_start:index_end,:])+ np.asarray(h5f[index_start+int(h5f.shape[0]/2):index_end+int(h5f.shape[0]/2),:]))/2.0
    snpEffects.append(snp_temp)


coor = pd.read_csv(args.coorFile,sep='\t',header=None,comment='#')
coor = coor.iloc[index_start:index_end,:]

#comptue expression effects
n_snps = coor.shape[0]
snpExpEffects = compute_effects(snpEffects, \
                                n_snps,\
                                models, maxshift=maxshift, nfeatures=args.nfeatures,
                                batchSize = args.batchSize, old_format = old_format)
#write output
snpExpEffects_df = coor
snpExpEffects_df=pd.concat([snpExpEffects_df.reset_index(),pd.DataFrame(snpExpEffects, columns = modelList.iloc[:,1])],axis=1,ignore_index =False)
snpExpEffects_df.to_csv(args.output, header = True)







