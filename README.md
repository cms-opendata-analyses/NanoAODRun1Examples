# NanoAODRun1Examples

Public usage examples for NanoAODRun1 are being provided here to illustrate how CMS Open Data in a NanoAOD or NanoAOD-like format can be analyzed to (re)produce a histogram of a physics observable -- specifically, the dimuon mass spectrum. These examples highlight the fact that analysis code for NanoAOD files is the same irrespective of data taking year or era, and can be analyzed with almost any version of ROOT (or ROOT-compatible) software (see e.g. ROOT User's guide, the python interface PyROOT, the RDATAframe framework, uproot, etc. ...) with any modern operating system.

More information on these examples can also be found in this [CERN TWiki](https://twiki.cern.ch/twiki/bin/view/CMSPublic/NanoAODRun1Examples), and this [CMS Open data workshop presentation](https://indico.cern.ch/event/1139022/contributions/4969980/attachments/2488290/4273320/NanoAODRun1ODWS.pdf).

CMS anticipates that about 50% of its publishable physics analyses can be performed using information stored in the NanoAOD format. For Run 1 Open Data, the NanoAOD format can be mimicked by processing AOD files through the [NanoAODRun1ProducerTool](https://opendata.cern.ch/record/12505). This tool is not intended to be pedagogical, but can be run using the CMSSW docker containers or virtual machines. The NanoAOD files needed for these analysis examples are available on the Open Data Portal as derived datasets.

Users interested in pedagogical code at AOD level are referred to the [POET framework](https://opendata.cern.ch/record/12501), which can be used to create NanoAOD-like files with content customized by the user. 

## Installation

The best way run these examples is to access the ROOT docker container described on the [Open Data Portal docker guide](https://opendata.cern.ch/docs/cms-guide-docker). You can also install [ROOT](https://root.cern.ch) locally on your computer. You will need a relatively recent ROOT version and the [XRootD package](https://xrootd.slac.stanford.edu/index.html) installed.

Clone this repository to access the example scripts:

```
$ git clone https://github.com/cms-opendata-analyses/NanoAODRun1Examples.git
$ cd NanoAODRun1Examples/
```

In the ROOT docker container, X11 graphics can be viewed using the built in VNC connection. In the container, execute:
```
$ start_vnc
Starting VNC server...
[1] 467
[2] 471
[3] 472
[4] 473
VNC connection points:
        VNC viewer address: 127.0.0.1::5901
        HTTP access: http://127.0.0.1:6080/vnc.html
```
Paste the HTTP access URL into your browser (or the VNC viewer address into a standalone VNC program) and select "connect". Provide the default password `cms.cern`. All images from the ROOT scripts in these examples will appear on this VNC screen. 

At the end of your session, before closing the docker container, stop the VNC connection by executing `stop_vnc`.

## 2010 Open Data example

Create a histogram of the dimuon mass spectrum in 2010 data by running:

```
$ cd dimuon_2010/
$ root -l -b -q MuHistos_eospublic.cxx++

Processing MuHistos_eospublic.cxx...
reading root://eospublic.cern.ch//eos/opendata/cms/derived-data/NanoAODRun1/01-Jul-22/Run2010B_Mu_merged.root
writing to MuHistos_Mu_eospublic.root
entries = 26718043
event nr 0
...periodic event number update...
```

Depending on network connection, this might take about 40 minutes. The script will produce a ROOT file containing several histograms. Reading the script will show you how to:

* Open ROOT files over the network
* Create histograms
* Extract a data tree from the input file
* Loop over the tree, accessing one event at a time
* Make event selections using branches in the tree
* Fill histograms with physics observables
* Write the histograms to an output file

To draw the mass spectrum:
```
$ start_vnc # only if not done already in this session
$ root -l MuHistos_Mu_eospublic.root
root [1] GM_mass_log->Draw()
root [2] c1->SetLogy(1)
root [3] .q
```

![dimuon histogram from 2010 in log-log scale](dimuon_2010/MuHistos_eospublic_mass.png)

This result is almost identical to a published result in [this paper](https://inspirehep.net/literature/1118729).

Another example script, called `MuHistos_publicchain.cxx++`, performs the same analysis, but shows how to read many individual files into a `TChain` using wildcards, rather than reading a single merged ROOT file.

## 2011 Open Data example

For 2011 data, NanoAODRun1 files can be used to demonstrate a basic analysis workflow using standard `TTree` selection methods or the [RDataFrame]() package. 

The TTree example script shows how the `TTree->Draw()` method can be used to apply event selections during histogram creation. It also shows how to create a plot canvas
with multiple histograms using varied plotting styles. 

To run the TTree example (this may take several hours, depending on network connection):
```
$ start_vnc # only if not done already in this session
$ cd dimuon_2011/
$ root -l Dimuon2011_eospublic.C

Processing Dimuon2011_eospublic.C...
high pt dimuon
all MuOnia
Quarkonium/Low pT dimuon only
Quarkonium and Upsilon
Quarkonium and B0
Quarkonium and Jpsi
Quarkonium and Jpsi/psiprime
```

You will see (briefly!) a plot similar to Figure 11 in [this Conference Report](https://inspirehep.net/literature/1292243)
![dimuon plot, no labels](dimuon_2011/Dimuon2011_eospublic.png)

The RDataFrame script produces the same plot in the end, but shows how to process a NanoAOD file through RDataFrame actions such as `Filter()`, `Define()`, and `Histo1D()`. 

To run the RDataFrame example (this should be much quicker, 10-20 minutes), first determine how many threads are accessible on your machine. If you wish to use fewer than 12 threads, edit the file `dimuon_2011/Dimuon2011_eospublic_RDF2.C` in a text editor and reduce the `nThreads` variable to a smaller number.
```
$ start_vnc # only if not done already in this session
$ cd dimuon_2011/
$ root -l Dimuon2011_eospublic_RDF2.C
Processing Dimuon2011_eospublic_RDF2.C...
Info in <TCanvas::Print>: file Dimuon2011_eospublic_RDF2.png has been created
Info in <TCanvas::Print>: pdf file Dimuon2011_eospublic_RDF2.pdf has been created
Run number: pass=24572802   all=35287778   -- eff=69.64 % cumulative eff=69.64 %
Dimuon threshold: pass=5277068    all=24572802   -- eff=21.48 % cumulative eff=14.95 %
Run number: pass=40550056   all=55915322   -- eff=72.52 % cumulative eff=72.52 %
Dimuon threshold and sample overlap: pass=40206351   all=40550056   -- eff=99.15 % cumulative eff=71.91 %
HLT       : pass=3418759    all=40206351   -- eff=8.50 % cumulative eff=6.11 %
Run number: pass=40550056   all=55915322   -- eff=72.52 % cumulative eff=72.52 %
Dimuon threshold & sample overlap: pass=40206351   all=40550056   -- eff=99.15 % cumulative eff=71.91 %
J/Psi threshold: pass=35226938   all=40206351   -- eff=87.62 % cumulative eff=63.00 %
HLT       : pass=23438007   all=35226938   -- eff=66.53 % cumulative eff=41.92 %
Run number: pass=40550056   all=55915322   -- eff=72.52 % cumulative eff=72.52 %
Dimuon threshold and sample overlap: pass=40206351   all=40550056   -- eff=99.15 % cumulative eff=71.91 %
Run number: pass=40550056   all=55915322   -- eff=72.52 % cumulative eff=72.52 %
Dimuon threshold and sample overlap: pass=40206351   all=40550056   -- eff=99.15 % cumulative eff=71.91 %
Run number: pass=40550056   all=55915322   -- eff=72.52 % cumulative eff=72.52 %
Dimuon threshold and sample overlap: pass=40206351   all=40550056   -- eff=99.15 % cumulative eff=71.91 %
Run number: pass=40550056   all=55915322   -- eff=72.52 % cumulative eff=72.52 %
Dimuon threshold and sample overlap: pass=40206351   all=40550056   -- eff=99.15 % cumulative eff=71.91 %
Run number: pass=40550056   all=55915322   -- eff=72.52 % cumulative eff=72.52 %
Dimuon threshold & sample overlap: pass=40206351   all=40550056   -- eff=99.15 % cumulative eff=71.91 %
HLT       : pass=10891153   all=40206351   -- eff=27.09 % cumulative eff=19.48 %
Elapsed time in seconds: 1061 sec
Elapsed time in minutes: 17 min
Number of threads: 12
```

You will see a plot similar to Figure 68 in [this paper](https://inspirehep.net/literature/1485699).
![dimuon plot, labels](dimuon_2011/Dimuon2011_eospublic_RDF2.png)

To close the ROOT session after viewing either plot, enter `.q` at the ROOT prompt. Remember to use `stop_vnc` to close the docker container graphics connection when you are finished with your session.

## 2012 Open Data example

NanoAOD-like outreach files have been produced for 2012 data in a very lightweight format that contains only the muon information. An outreach example analysis using these files can be found at [this Open Data Portal record](https://opendata.cern.ch/record/12342). However, the same type of analysis can be performed using full NanoAODRun1 (or newer NanoAOD) files. 

2012 examples using NanoAODRun1 files exist in both ROOT C++ and PyROOT forms. The scripts are functionally identical and both use the RDataFrame method for quicker speeds, but show the differences between using RDataFrame in C++ versus Python. 

To run the C++ script:
```
$ start_vnc # only if not done already in this session
$ cd dimuon_2012/
$ root -l dimuonSpectrum2012_eospublic.C
root [0]
Processing dimuonSpectrum2012_eospublic.C...
Info in <TCanvas::Print>: pdf file dimuonSpectrum2012_C_eospublic.pdf has been created
Events with two or more muons: pass=59571658   all=61540413   -- eff=96.80 % cumulative eff=96.80 %
Muons with opposite charge: pass=37597805   all=59571658   -- eff=63.11 % cumulative eff=61.09 %
root [1] .q
```

To run the Python script:
```
$ start_vnc # only if not done already in this session
$ cd dimuon_2012/
$ python dimuonSpectrum2012_eospublic.py
Info in <TCanvas::Print>: pdf file dimuonSpectrum2012_py_eospublic.pdf has been created
Events with at least two muons: pass=59571658   all=61540413   -- eff=96.80 % cumulative eff=96.80 %
Muons with opposite charge: pass=37597805   all=59571658   -- eff=63.11 % cumulative eff=61.09 %
```

You will see a plot like this:
![2012 dimuon plot](dimuon_2012/dimuonSpectrum2012.png)

Remember to use `stop_vnc` to close the docker container graphics connection when you are finished with your session.

## Downloading files locally

All of these examples use the XRootD protocol to stream the data files over your network connection. If you prefer to download the files locally (you'll need some disk space!)
this can be done using `xrdcp`. In the ROOT container, a command like this will download the stated file (from the 2010 example) locally:

```
$ xrdcp root://eospublic.cern.ch//eos/opendata/cms/derived-data/NanoAODRun1/01-Jul-22/MuHistos_Mu_eospublic.root .
```

This download may be slow, depending on your network connection and the size of the files, but only needs to be done once. 
To use the local files, the `TFile::Open` commands in the example scripts should be edited to replace the entire `root://eospublic.cern.ch/..../FILE.root` file name with the local file path. 
