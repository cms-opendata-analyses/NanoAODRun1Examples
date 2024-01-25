# NanoAODRun1Examples

Public usage examples for NanoAODRun1 are being provided here to illustrate how CMS Open Data in a NanoAOD or NanoAOD-like format can be analyzed to (re)produce a histogram of a physics observable -- specifically, the dimuon mass spectrum. These examples highlight the fact that analysis code for NanoAOD files is the same irrespective of data taking year or era, and can be analyzed with almost any version of ROOT (or ROOT-compatible) software (see e.g. ROOT User's guide, the python interface PyROOT, the RDATAframe framework, uproot, etc. ...) with any modern operating system.

More information on these examples can also be found in this CERN TWiki, and this CMS Open data workshop presentation.

CMS anticipates that about 50% of its publishable physics analyses can be performed using information stored in the NanoAOD format. For Run 1 Open Data, the NanoAOD format can be mimicked by processing AOD files through the NanoAODRun1ProducerTool. This tool is not intended to be pedagogical, but can be run using the CMSSW docker containers or virtual machines. The NanoAOD files needed for these analysis examples are available on the Open Data Portal as derived datasets.

Users interested in pedagogical code at AOD level are referred to the POET framework, which can be used to create NanoAOD-like files with content customized by the user. 

## Installation

The best way run these examples is to access the ROOT docker container described on the Open Data Portal docker guide. You can also install ROOT locally on your computer. You will need a relatively recent ROOT version and the XRootD package installed.

Clone this repository to access the example scripts:

```
$ git clone https://github.com/cms-opendata-analyses/NanoAODRun1Examples.git
$ cd NanoAODRun1Examples/
```

## 2010 Open Data example

Create a histogram of the dimuon mass spectrum in 2010 data by running:

```
$ cd dimuon_2010/
$ root -l MuHistos_eospublic.cxx++
```

Depending on network connection, this might take about 40 minutes. The script will produce a ROOT file containing several histograms. To draw the mass spectrum:

```
$ start_vnc # for viewing an X11 window from docker
$ root -l MuHistos_Mu_eospublic.root

root [1] GM_mass_log->Draw()
```

![dimuon histogram from 2010 in log-log scale](MuHistos_eospublic.png)




