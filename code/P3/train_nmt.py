import numpy
import os
import sys

from nmt import train

def main(job_id, params):
    print params
    validerr = train(saveto=params['model'][0],
                     reload_=params['reload'][0],
                     dim_word=params['dim_word'][0],
                     dim_enc=params['dim_enc'], # multi layer
                     dim_dec=params['dim_dec'][0], 
                     dim_coverage=params['dim_coverage'][0],
                     down_sample=params['down_sample'],
                     dim_target=params['dim_target'][0],
                     dim_feature=params['dim_feature'][0],
                     decay_c=params['decay-c'][0],
                     clip_c=params['clip-c'][0],
                     lrate=params['learning-rate'][0],
                     gamma=params['alphas-gamma'][0],
                     optimizer=params['optimizer'][0], 
                     patience=15,
                     maxlen=params['maxlen'][0],
                     batch_size=8,
                     valid_batch_size=8,
                     validFreq=-1,
                     dispFreq=100,
                     saveFreq=-1,
                     sampleFreq=-1,
          datasets=['../data/online-train.pkl',
                    '../data/train_caption.txt',
                    '../data/align-online-train.pkl'],
          valid_datasets=['../data/online-test.pkl',
                    '../data/test_caption.txt'],
          dictionaries=['../data/dictionary.txt'],
                         use_dropout=params['use-dropout'][0])
    return validerr

if __name__ == '__main__':
    
    maxlen=[400]
    modelDir='models/'
    dim_word=[256]
    dim_enc=[250,250,250,250] # they are bidirectional
    dim_dec=[256]
    dim_coverage=[121]
    down_sample=[0,0,1,1]

        
    main(0, {
        'model': [modelDir+'attention_maxlen'+str(maxlen)+'_dimWord'+str(dim_word[0])+'_dim'+str(dim_dec[0])+'.npz'],
        'dim_word': dim_word,
        'dim_enc': dim_enc,
        'dim_dec': dim_dec,
        'dim_coverage': dim_coverage,
        'down_sample':down_sample,
        'dim_target': [111], 
        'dim_feature': [8],
        'optimizer': ['adadelta'],
        'decay-c': [0.], 
        'clip-c': [1000.], 
        'use-dropout': [False],
        'learning-rate': [1e-8],
        'alphas-gamma': [0.1],
        'maxlen': maxlen,
        'reload': [False]})



