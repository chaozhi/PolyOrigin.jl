# CHANGELOG


# TODO list 

* [] error if eg. GT = "1/1/1" (line198, likelihiood.jl)
* [] error/war if poid is not specfied
* [] check results for offspring with all data being missing
* [] The marker error plot does not have the chromosomes in order. 
* [] tutorial for using polyorigin in an updating way
* [] docs
* [] test for 6x, accuracy of ancestry inference. 

* [x] test byneighbor = true
* [x] reset isrefinedis = nothing (true only if size>=50)
* [x] reset tukeyfence = nothing (3 if size>=50 and Inf otherwise)
* [x] debug error in skeloton distance
* [x] allow byparent = 2
* [x] simulate data with a multi-generation pedigree
* [x] polyreconstruct viterbi
* [x] link offspring to funders for FB or viterbi results
* [x] plotviterbipath
* [x] remove byparent, using only byparent = true
* [x] reject parentcorrect if not increasing likelihood; set parentcorrect=false in polyphase
* [x] polyphase calculate logl accounting for prior prob of parent genotypes
* [x] plotcondprob for midparent
