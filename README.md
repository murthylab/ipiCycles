# supplemental data and code for Stern et al. 2017

## data
All raw data are in `data/` and contain pulse times and ipis. The files starting with `dat/Stern2014*` contain different annotations of the song recordings (Canton S, period mutants, D. simulans) first published in Stern et al. 2014 (see David Stern's [website][1] for the raw recordings):
- `dat/Stern2014_FSSStern2014.mat` - automatically segmented using the original version of [FlySongSegmenter][2] and taken from the [data supplement][1] to [Arthur et al.,2013][3].
- `dat/Stern2014_FSSStern2014.mat` - automatically segmented using a modified version of the [FlySongSegmenter from the Murthy lab][4] with parameters used in [Coen et al. (2014)][5]
- `dat/Stern2014_KyriacouManual2017.mat` - manual annotations of parts of the recordings (Canton S and perL) downloaded from the data supplement to [Kyriacou et al. (2017)][6]. `dat/Stern2014_KyriacouManual2017_pulses.mat` contains the pulse waveforms for the manually annotated pulses.

`dat/Stern2017_FSSCoen2014` contains segmentations of new [recordings][1] of Canton S and perL flies segmented using the [Murthy lab segmenter][4]

## code
`Fig*` folders contain code and data to reproduce the figures in the paper - see comments in the respective scripts. `code` contains general source code used by the scripts. `spectra` contains code and data for pre-calculating the IPI spectra for all data sets. 


[1]: https://www.janelia.org/lab/stern-lab/tools-reagents-data "Stern lab data"
[2]: https://github.com/FlyCourtship/FlySongSegmenter "original fly song segmenter"
[3]: [1]: https://bmcbiol.biomedcentral.com/articles/10.1186/1741-7007-11-11 "Multi-channel acoustic recording and automated analysis of Drosophila courtship songs"
[4]: https://github.com/murthylab/songSegmenter "Murthy lab song segmenter"
[5]: https://www.nature.com/nature/journal/v507/n7491/full/nature13131.html "Dynamic sensory cues shape song structure in _Drosophila_"
[6]: http://www.pnas.org/content/114/8/1970.abstract "Failure to reproduce period-dependent song cycles in Drosophila is due to poor automated pulse-detection and low-intensity courtship"
