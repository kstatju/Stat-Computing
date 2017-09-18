####################################################################################
# Problem 3 



likelytest <-  function(sample = NULL, dist = "Poisson", dist.with.pram = 1, null.hypo = 1, sample_size = 1000, 
                        sim_p.value = FALSE, replicat.times = 10000){
    
    if (sim_p.value){
      if (!is.null(sample)) dist.with.pram = mean(sample)
      
      likelyhoodlist = replicate(replicat.times, likelytest(sample = NULL, dist = dist, dist.with.pram = null.hypo, 
                                                            null.hypo = null.hypo, sample_size = sample_size))
    }
    if (is.null(sample)){
        if (toupper(dist) == "POISSON") { rand = rpois; pand = dpois    }
        rand_var = rand(n = sample_size, lambda = dist.with.pram)
    } else {rand_var = sample}
  
    alter_hypo = mean(rand_var)
    
    likelyhood = prod(pand(x = rand_var, lambda = null.hypo))/prod(pand(rand_var, lambda = alter_hypo))
    
  
  if (sim_p.value) {
        quantile_val = quantile(likelyhoodlist, probs = c(.95))
        return(c(test.statistics = likelyhood, Quantile.upper = quantile_val, 
                 p.value = 1 - sum(likelyhoodlist > likelyhood)/length(likelyhoodlist)))

    
    }else {
        return(test.statistics = likelyhood)
    }
}

kk = likelytest(dist.with.pram = 2, null.hypo = 2, sample_size = 10, sim_p.value = TRUE)
hist(kk)
kk
