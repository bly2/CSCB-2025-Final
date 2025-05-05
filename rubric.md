### Task 1: Create a method that infers CNAs from scRNA-seq data (20pts – 447 sections 10pts – 647 section)
•	Working Code in python that takes an input and minimally provides outputs in the required format
•	Output has to have at least the following a list of CNAs where each CNA is defined by genomic region (chromosome, start, stop,) and some indication of the type of CNA (e.g. 'gain' | loss') an assignment of CNAs to cells, or to groups of cells. (It's Ok to add to AnnData)
•	Code can be downloaded from GitHub and run. Provide a ReadMe file on how to do this
•	Method novelty – You should comment on the specific aspect of your method that is novel

### Task 2A: Assessment (10pts – 447 section; 10pts – 647 section)
•	Analyze your method's performance in terms of standard metrics such as accuracy, precision or area under precision recall curves on the test data. 
•	Provide a visualization/table of this
•	Compare the performance of your method to at least one other existing method. Provide a visualization/table for this
•	Re-run your evaluation for different parameters
•	If your method has a confidence score, evaluate this as well.  (optional - depends on your method)

### Task 2B: Augment assessment with better gold standard data (0pts – 447 section; 10pts – 647 section)
•	Obtain data that has CNAs of small, medium and large sizes. Similarly, ensure that the data has a range of CNAs present at low, medium and high frequencies
•	Inclusion of the source/citation of the data, whether from literature or from simulations. Provide a few sentences of how CNAs were identified in this 'gold standard data'
•	Include pre-processing of the data in your GitHub repo
•	Repeat the evaluations done in Task 2A

### Task 3: Measure CNA in PSCs (10pts – 447 section; 10pts – 647 section)
•	Code showing pre-processing steps of at least 3 PSC datasets
•	Tables/Visualizations comparing known CNAs and potentially new ones

### Task 4: Predict CNA impact (Extra Credit) (10pts – 447 section; 10pts – 647 section)
•	Brief description of the approach for predicting functional impact. For example, through GSEA, perturbations, etc.
•	Code assessing functional impact
•	Table/Visualization summarizing the results

### Paper (a total of 40 pts)						
•	Introduction: clearly stated CNA problem, background, motivation (4pts)
•	Methodology: Precisely described your developed tool and justified the design choices, partly by surveying existing methods (7pts)
•	Benchmarking: assess your method on data (simulated or real data) and show performance metric. Benchmark against at least 1 other method. (7pts)
•	Analysis: apply your method onto 3 or more iPSC data and discuss the previously reported PSC CNAs and new ones that you discovered. (8pts)
•	Discussion: depth of discussion. Future directions. (4pts)
•	Figures: 2~5 large figures with multiple panels. Clearly written captions (5pts)
•	Organization: 2000~5000 words. Clarity of writing. (5pts)

### Presentation (a total of 20 pts)
•	Content: clearly explained the problem, methodology, benchmark, analysis, and conclusion (6pts)
•	Organization: stay within time limit, logically structured (3pts)
•	Visual: effective usage of figures/tables to support the narrative (3pts)
•	Delivery: confident, clear speech; well-prepared (2pts)
•	Q&A: thoughtful responds (6pts)

### Assessment of Teammate (a total of 6 pts)
Please score each team member from 1 (superlative) to 10 (unacceptable): 
•	Contribution to planning (1pt)
•	fulfill assigned task (1pt)
•	timeliness (1pt)
•	responsiveness to communication (1pt)
•	pro-activeness (1pt)
•	1-3 sentences of each teammate's role (1pt)

### Assess Another Team (a total of 14 pts)
Score your assigned team's code & paper. 
Paper: 
•	did it address task 1 (1pt)
•	did it address task 2 (1pt)
•	did it address task 3 (1pt)
•	clarity of writing (1pt)
•	depth of discussion (1pt)
•	introduction (1pt)
•	figures (1pt)
•	1 paragraph of paper summary (2pts)
Code: 
•	ease of installation (1pt)
•	clarity of documentation (1pt)
•	ease of running based on README (1pt)
•	experience when applying to a new dataset (1pt)
•	reproducibility (1pt)


