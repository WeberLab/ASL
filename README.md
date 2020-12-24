# ASL

This was mostly written by Mike Jarrett (mike.jarrett@ubc.ca) in 2016-2018, with contributions from Alex Weber and Esther Lin.

Alex Weber modified the code for a 3d m0

## Dependencies

You must have python3 install to run this.

* numpy
* nibabel
* nipy

## Using this package

Create a directory where you keep your git repositories. For this, we'll use `~/repos/`.

```
mkdir ~/repos/
cd ~/repos/
git clone git@github.com/WeberLab/ASL.git
```

### ASL

To create a CBF map from an ASL scan:

`~/repos/ASL/asl.py [-p PLD] [-s SliceDelay] [-l LabelDuration] m0.PAR pCASL.PAR`

Important: PLD, SliceDelay and LabelDuration have default values of 
    PLD=1.60
    SliceDelay=0.039
    LabelDuration=1.65
These were current for the 3T protocol as of late 2017 but PLEASE CHECK THEM FOR YOUR OWN STUDY.
