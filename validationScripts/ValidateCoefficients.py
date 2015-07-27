import pandas as pd
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

def main(args):
    trueCoeffs = pd.read_table(args.trueCoeffs, sep=',')

    def joinNames(X, suff=''):
        print("{}-{}".format(X['IID1'+suff], X['IID2'+suff]))
        return "{}-{}".format(X['IID1'+suff], X['IID2'+suff])

    trueCoeffs[['IID1', 'IID2']] = trueCoeffs[['IID1', 'IID2']].astype(str)
    trueCoeffs['key'] = trueCoeffs.apply(joinNames, axis=1)
    trueCoeffs.set_index('key', inplace=True)

    estCoeffs = pd.read_table(args.estCoeffs, sep='\t',
                              skiprows=1,
                              names=['IID1Est', 'IID2Est',\
                                     'k0_hat', 'k1_hat', 'k2_hat',\
                                     'Estimated Kinship Coefficient',\
                                     'nbSNP'])

    estCoeffs[['IID1Est', 'IID2Est']] = estCoeffs[['IID1Est', 'IID2Est']].astype(str)
    def joinNamesEst(X):
        return joinNames(X, suff='Est')

    estCoeffs['key'] = estCoeffs.apply(joinNamesEst, axis=1)
    estCoeffs.set_index('key', inplace=True)
    jointCoeffs = trueCoeffs.join(estCoeffs)

    sns.set_context('poster')
    sns.set_style('whitegrid')
    plt.scatter(jointCoeffs['Expected Kinship Coefficient'].values,
                jointCoeffs['Estimated Kinship Coefficient'].values)
    plt.xlabel("Expected Kinship Coefficients")
    plt.ylabel("Estimated Kinship Coefficients")
    sns.despine()
    plt.show()
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Validator")
    parser.add_argument('--trueCoeffs', help="file with true coefficients",
                        type=argparse.FileType('r'))
    parser.add_argument('--estCoeffs', help="file with estimated coefficients",
                        type=argparse.FileType('r'))

    main(parser.parse_args())
