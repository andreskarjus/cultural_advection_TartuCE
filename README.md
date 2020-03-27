# "Selection above the baseline: quantifying the advection effect in four domains of cumulative culture"

- Scripts to replicate the analyses presented at the [Applications in Cultural Evolution](https://cultevol.ut.ee/) conference, June 5-8, Tartu, Estonia.
- All the scripts are in R. You'll need a bunch of packages; sourcing the functions file will attempt to install the missing ones; alternatively, run this before running any scripts: 
`install.packages(c("visNetwork", "igraph", "XML", "magrittr", "plotly", "crosstalk", "parallel", "fastmatch", "text2vec", "grr", "Matrix", "compiler"))`
- See the slides https://andreskarjus.github.io/cultevol_tartu_slides for more information. 
- Check the scripts for the places to download the datasets (or the handy unofficial dumps).
- Code to apply the advection model to language data (the Corpus of Historical American English in particular) is available here: https://github.com/andreskarjus/topical_cultural_advection_model

- If you make use of these scripts, please cite the following paper which defines and describes the model: _Andres Karjus, Richard A. Blythe, Simon Kirby, Kenny Smith, 2018. Quantifying the dynamics of topical fluctuations in language. arXiv preprint. https://arxiv.org/abs/1806.00699_ --> UPDATE the paper is now published in Language Dynamics and Change: https://doi.org/10.1163/22105832-01001200, please cite that instead; the preprint remains freely available as well.



- Finally, if you find the customized plotting scripts useful and adopt them for your own work, besides the obvious package citations, I'd appreciate a mention of the following repo: Karjus, Andres (2018). aRt of the Figure. GitHub repository, https://github.com/andreskarjus/artofthefigure. bib:
```
@misc{karjus_artofthefigure_2018, 
  author = {Karjus, Andres}, 
  title = {aRt of the Figure}, 
  year = {2018}, 
  publisher = {GitHub}, 
  journal = {GitHub repository}, 
  howpublished = {\url{https://github.com/andreskarjus/artofthefigure}},
  DOI = {10.5281/zenodo.1213335}
} 
```

- Feel free to get in touch if something does not work.
