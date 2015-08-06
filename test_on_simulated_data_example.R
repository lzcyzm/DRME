
res <- simulateData(n_Sites=20000,replicate=3,min_expression=1,max_expression=4,
                    dif_express=1/2,per_me=1/2,dif_me=1,lib_s=1/5,var=0.2)


res <- DMEseq(res[[1]],res[[2]],res[[3]],res[[4]])