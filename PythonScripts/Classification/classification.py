import argparse
import numpy as np
import pandas as pd
from sklearn.model_selection import cross_val_score
#from svm_l0 import SVM_L0
from sklearn.svm import SVC, LinearSVC

def ReadMDFile(file):
    with open(file) as fp:
        lines = [line.rstrip('\n') for line in fp]
    height = len(lines)
    width = len(lines[0].split(' ')) - 1
    print(height, width)
    md_mat = np.zeros((height, width))

    for y in xrange(height):
        elems = lines[y].split(' ')
        for x in xrange(width):
            md_mat[y, x] = float(elems[x])

    return md_mat


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('-md_file', dest='md_file', type=str,
                        help='path to the md matrix')
    parser.add_argument('-labels', dest='labels', type=str,
                        help='path to the labels file')
    args = parser.parse_args()
    md = ReadMDFile(args.md_file)
    print(md)
    labels = ReadMDFile(args.labels).ravel()
    print(labels)
    print(labels.shape)
    clf = LinearSVC(C=10000.0)#SVM_L0(verbose=1, feature_selection=True)
    splits = min((labels==1).sum(), (labels==-1).sum())
    clf1 = LinearSVC(C=10000.0)
    clf1.fit(md, labels)
    print(clf1.score(md, labels))
    scores = cross_val_score(clf, md, labels, cv=splits, scoring='balanced_accuracy')
    print(scores)
    print(list(scores))
