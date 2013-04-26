jam-pipeline
============

This is a development project of the Putnam Lab at Rice University.


Quick start
=================

How to download and run the tests.  Two test datasets are available, both derived from the Limulus polyphemus JAM project:  
https://dl.dropboxusercontent.com/u/8774493/JAMtestSequenceQuick.tgz and https://dl.dropboxusercontent.com/u/8774493/JAMtestSequence.tgz
	
  - cd ~; git clone https://github.com/putnamlab/jam-pipeline.git
  - Download the test data from https://dl.dropboxusercontent.com/u/8774493/JAMtestSequenceQuick.tgz 
  - Uncompress the archive:  cd ~;  gzcat JAMtestSequenceQuick.tgz | tar -xvf -
  - export JAM_ROOT=~/JAMtestSequenceQuick
  - cd ~/jam-pipeline
  - export PATH=~/jam-pipeline/scripts:~/jam-pipeline/bin:$PATH
  - make all
  - make test
  - look at the files and directories that were created in $JAM_ROOT/JAMtestSequenceQuick/Limulus_testpolyphemus/



Licence
=========


JAM Pipeline is available at https://github.com/putnamlab/jam-pipeline and is realeased under the [Creative Commons Attribution 3.0 Unported License](http://creativecommons.org/licenses/by/3.0/deed.en_US).

<a rel="license" href="http://creativecommons.org/licenses/by/3.0/deed.en_US">
<img alt="Creative Commons License" style="border-width:0" src="http://i.creativecommons.org/l/by/3.0/88x31.png" />
</a>


Credits 
=========

by Paul Havlak and Nicholas Putnam
