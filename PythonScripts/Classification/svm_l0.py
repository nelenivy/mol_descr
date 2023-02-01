# ========================================================
# SVM with zero norm from Alexey Shestov <github: nelenivy>
#    Realization of Jason Weston "Use of the Zero-Norm with Linear Models and Kernel Methods"
#    Approximation of the zero-norm Minimization (AROM) - add multiplicative update
#    Use l2-SVM => l2-AROM
#    Regularization of linear inseparable case:
#       Make linear separation by hands - add new columns to kernel matrix
#
#    sci-kit learn trainer
#    
#    SVM_L0 - class for binary classification
# ========================================================

import sys
import numpy as np

from sklearn.base import BaseEstimator, ClassifierMixin, clone
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC, LinearSVC


class SVM_L0(BaseEstimator, ClassifierMixin):
    def __init__(self, scale=True, C=1.0, eps=1.0e-3, \
                 estimator=LinearSVC(C=100000.0), verbose=0, \
                 feature_selection=False, predict_estimator=SVC(C=1.0, kernel='linear'), \
                 n_iter=1000):
        """
        feature_selection:
            True: use predict_estimator to predict data
            False: use 'z' coefs to predict data
        scale:
            don't supported, always acts like scale=True: uses StandardScaler
        """

        self.verbose = verbose

        # use algorithm for feature selection or for training
        self.feature_selection = feature_selection
        self.predict_estimator = predict_estimator
        self.trainer = None

        self.scale = scale
        self.scaler = StandardScaler()
        # parameter for regularization
        self.C = C

        # stop criterion parameters
        self.eps = eps
        self.n_iter = n_iter

        # estimator
        self.estimator = estimator

        # parameters for fit
        self.z = None
        self.inds = None
        self.intercept_ = None

        return

    def prepareGraphics(self):
        """
        """
        from matplotlib import pyplot as plt

        # initialize
        plt.ion()
        self.fig = plt.figure()
        self.ax1 = self.fig.add_subplot(121)
        self.ax2 = self.fig.add_subplot(122)

        return

    def drawZ(self, z_full, z_old_full, ymax=2.0):
        """
        """

        self.ax1.clear()
        self.ax2.clear()
        # self.ax1.set_ylim(0,ymax)
        self.ax1.bar(range(len(z_full)), z_full, color='blue')
        self.ax2.set_ylim(0, 0.1)
        self.ax2.bar(range(len(z_full)), abs(z_full - z_old_full))
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

        return

    def fit(self, X, y):
        """
        Make SVML0 without matrix update (without removing unused features)
        """

        if self.verbose >= 10:
            self.prepareGraphics()

        y = np.array(y)

        # y must have +1,-1 labels, check
        assert np.all((y == 1) + (y == -1))
        self.non_zero_cols = (X > 0.0).sum(axis=0) > 0
        X = X[:, self.non_zero_cols]
        # copy X, scale it , add regularization cols
        X = self.scaler.fit_transform(X)

        # regularization matrix
        _reg_mat = 1.0 / self.C * np.ones((X.shape[0], X.shape[0]))

        # add regularization mat to data matrix
        X = np.concatenate((X, _reg_mat), axis=1)

        # selected columns
        # sel_inds = np.array(range(X.shape[1]))
        # sel_inds_prev = np.array(sel_inds)
        z = np.ones((1, X.shape[1]))
        z_prev = np.zeros((1, X.shape[1]))

        inds = np.array(range(X.shape[1]))

        trainer = None

        _counter = 0
        while not (np.all(abs(z - z_prev) < self.eps * np.max(z)) \
                   or _counter > self.n_iter):
            print _counter
            # get subset from previous subset
            subset = abs(z[0]) > self.eps * np.max(z)
            z = z[:, subset]
            inds = inds[subset]

            # prepare matrix mat
            mat = X[:, inds] * z

            # calculate new z
            trainer = clone(self.estimator)
            # set trainer to be weighted
            trainer.set_params(class_weight={1: float(np.sum(y == -1)) / float(np.sum(y == 1))})
            trainer.fit(mat, y)

            z_prev = z
            z = z * trainer.coef_

            if self.verbose > 0:
                sys.stdout.write("\r" + str(_counter) + " | " + \
                                 str(np.max(abs(z - z_prev))) + " | " + \
                                 str(np.max(z)) + " | " + str(len(inds)) + " | ")
                sys.stdout.flush()

                if self.verbose >= 10:
                    # draw z full
                    _z_full = np.zeros(X.shape[1])
                    _z_full[inds] = z[0]
                    _z_old_full = np.zeros(X.shape[1])
                    _z_old_full[inds] = z_prev[0]
                    self.drawZ(_z_full, _z_old_full)

            _counter += 1

            # remove regularizator inds
        subset = inds < (X.shape[1] - X.shape[0])
        inds = inds[subset]
        z = z[:, subset]

        self.inds = inds
        self.z = z
        self.intercept_ = trainer.intercept_

        if self.feature_selection:
            self.trainer = clone(self.predict_estimator)
            self.trainer.set_params(class_weight={1: float(np.sum(y == -1)) / float(np.sum(y == 1))})
            mat = X[:, self.inds]
            self.trainer.fit(mat, y)
            print self.trainer.score(mat, y)
        # self.z = z[:,:(X.shape[1]-X.shape[0])]

        if self.verbose > 0:
            sys.stdout.write("\n")

        return

    def predict(self, X):
        """
        """
        X = X[:, self.non_zero_cols]
        X = self.scaler.transform(X)

        if not self.feature_selection:
            res = np.matmul(X[:, self.inds], self.z.T) + self.intercept_
            prediction = (res > 0).astype(int) - (res <= 0).astype(int)
        else:
            mat = X[:, self.inds]
            prediction = self.trainer.predict(mat)

        return prediction
