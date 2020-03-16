# AMATH482 - Computational Methods in Data Analysis
Assignments for Computational Methods in Data Analysis class. Cover a range of topics including the Fast Fourier Tranform, Singular Value Decomposition, Prinicipal Component Analysis and some Machine Learning for applications in image/signal processing and classifications.

# Important Information/How to Use
- Many of the assignments require additional files, which have not been uploaded to this repository.
- Code for the first four assignments has been written in MATLAB. Code for the fifth assignment was written using Python.
- Each assignment folder contains the code and a writeup with theoretical background, algorithm implementation, and summary and     conclusion information.

# What does each assignment do?
## Homework 1 - Fast Fourier Transforms
Given a noisy dataset representing an ultrasound of a marble at twenty different points in time, I used Fast Fourier Transforms to de-noise the data in order to find the marble's postion at the 20th data point.
## Homework 2 - Creating and filtering Spectrograms
Using Fast Fourier Transforms and filters (ex. Gabor Filter), we created spectrograms to analyze two recordings of "Mary Had a Little Lamb" played on a recorder and a piano. The goal was to re-create the score for the piece on both instruments after filtering out overtones and other unnecessary frequencies.
## Homework 3 - Singular Value Decomposition
Twelve different videos of a paint can moving up and down were provided. There were four cases, with three cameras recording per case at different angles in the room. The first two cameras recorded motion vertically, while the third camera recorded the motion horizontally. The test cases were as follows:
- Test 1: Ideal- no shaky camera
- Test 2: Shaky camera
- Test 3: Additional horizontal motion to the paint can
- Test 4: Additional horizontal and rotational motion to the paint can
Our goal was to obtain the position of the paint can in each time frame for every video and perform Singular Value Decompositon on all of these tests to determine what it told us about the motion of the paintcan.
## Homework 4 - Machine Learning (Classification)
Using Principal Comoponent Analysis and Linear Discriminant Analysis, we built and trained a classifier to work with these three tests:
- Test 1: Identify the song's artist from one of three artists across different genres.
- Test 2: Identify the song's artist from one of three artists in the same genre.
- Test 3: Identify the song's genre from one of three genres.
## Homework 5 - Machine Learning (Neural Networks)
Using the Fashion-MNIST dataset, we optimized a Fully Connected Dense Neural Network and a Convolutional Neural Network by exploring how changing certain parameters influenced the accuracy of the classification.
These parameters inclued zero-padding, filter size, activation function, and the learning rate. This was done using the Tensflow package in Python.
