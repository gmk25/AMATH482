# AMATH482 - Computational Methods in Data Analysis
Assignments for Computational Methods in Data Analysis class. Cover a range of topics including Fast Fourier Tranforms, Singular Value Decomposition, Prinicipal Component Analysis and some Machine Learning for applications in image/signal processing and classifications.
Each assignment folder contains the code and a writeup with theoretical background, algorithm implementation, and summary and conclusion information.

Note- Many of the assignments require additional files, which have not been uploaded to this repository.
# Homework 1 - Fast Fourier Transforms
Given a noisy dataset representing an ultrasound of a marble at twenty different points in time, we used Fast Fourier Transforms to de-noise the data in order to find the marble's postion at the 20th data point.
# Homework 2 - Creating and filtering Spectrograms
Using Fast Fourier Transforms and filters (ex. Gabor Filter), we created spectrograms to analyze two recordings of "Mary Had a Little Lamb" played on a recorder and a piano in order to re-create the score for the piece.
# Homework 3 - Singular Value Decomposition
Twelve different videos of a paint can moving through space were provided. There were four groups, with three videos in each group - Camera 1 and 2 recorded vertical motion, while Cmaera 3 recorded the motion horizontally. While all of the tests recorded the same paintcan, some tests included different motions for the paintcan and introduced shake to the camera. Our goal was to obtain the position of the paint can for each video and perform Singular Value Decompositon on each of these tests to determine what it told us about the motion of the paintcan.
# Homework 4 - Machine Learning (Classification)
We built and trained a classifier to work with these three tests:
Test 1: Identify the song's artist from one of three artists across different genres.
Test 2: Identify the song's artist from one of three artists in the same genre.
Test 3: Identify the song's genre from one of three genres.
This assignment also utilized Principal Component Analysis and Linear Discriminant Analysis.
# Homework 5 - Machine Learning (Neural Networks)
Using the Fashion-MNIST dataset, we optimized a Fully Connected Dense Neural Network and a Convolutional Neural Network by exploring how changing certain parameters influenced the accuracy of the classification.
These parameters inclued zero-padding, filter size, activation function, and learning rate. This was done using the Tensflow package in Python.
