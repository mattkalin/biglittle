# This algorithm was created by Matt Kalin in April 2020 
# It was used to create big/little pairings for Sigma Alpha Mu fraternity at 
# the University of Maryland, College Park during the Spring 2021 semester 
# Contact: kalinmatt4@gmail.com 
# git repository: https://github.com/mattkalin/biglittle


{
  ProbVector = function(x){
    return(x / sum(x))
  }
  AddCols = function(df, new.cols, default.val = NA){
    df[, new.cols] = default.val
    return(df)
  }
  MakePairs = function(people, groups, loss.fn, sample.size = 50,
                       n.sims = 5, identical.groups = FALSE, 
                       default.groups = NULL, alpha = 0.05, group.names = NA){

    # people is a character vector of everyone's names 
    # loss.fn is the loss function (high output is bad, low output is good)

    start.time <<- Sys.time()
    library(dplyr)
    n.groups = length(groups)
    if(n.groups == 1){
      n.groups = groups
      groups = rep(ceiling(length(people)/n.groups), n.groups)
    } 
    if(sum(groups) < length(people)){
      if(max(groups) == 1){
        stop("Too many people for this many groups. Some brothers need to take twins.")
      } else {
        stop("Too many people for this many groups. More brothers need to take twins.")
      }
      
    }
    if(isTRUE(is.na(group.names))){
      col.names = paste("Group", 1:n.groups)
    } else {
      col.names = group.names
    }
    prob.df = data.frame("Person" = people) %>% 
      AddCols(col.names, 1/n.groups)
    if(!is.null(default.groups)){
      for (p in names(default.groups)) {
        p.grp = default.groups[[p]]
        p.index = match(p, people)
        prob.df[p.index, col.names] = 0
        prob.df[p.index, col.names[p.grp]] = 1 / length(p.grp)
      }
    }
    default.prob.df <<- prob.df
    answers = vector('list', n.sims)
    best.pairs = vector('list', 2)
    best.pairs$Loss = Inf
    for (k in 1:n.sims) {
      sim.num <<- k
      print(paste("Simulation", k))
      answers[[k]]$Pairs = 
        MakePairsAux(groups, people, loss.fn, prob.df, n.groups, sample.size, 
                     1, identical.groups, alpha, shuffled.ppl = 
                       ifelse(rep(identical.groups, length(people)), 
                              sample(people), people), col.names = col.names)
      if(!is.null(answers[[k]]$Pairs$Sim.data)){
        return(answers[[k]]$Pairs)
      }
      answers[[k]]$Loss = loss.fn(answers[[k]]$Pairs, groups, n.groups, col.names) 
      print(paste("Loss for this simulation is", signif(answers[[k]]$Loss, 6)))
      if(answers[[k]]$Loss < best.pairs$Loss){
        best.pairs$Loss = answers[[k]]$Loss
        best.pairs$Pairs = answers[[k]]$Pairs
      }
    }
    min.loss = min(all.pairs.data$Loss) 
    if(min.loss < best.pairs$Loss){
      print("There was a pairing with better loss than any simulated answer")
      min.index = match(min.loss, all.pairs.data$Loss)
      best.row = all.pairs.data[min.index, ]
      new.pairs = data.frame("Person" = people) %>% 
        AddCols(col.names, 0)
      for (i in 1:length(people)) {
        new.pairs[i, col.names[best.row[, i]]] = 1
      }
      best.pairs$Pairs = new.pairs
      best.pairs$Loss = min.loss
      # don't need to change best loss because only pairs are returned
    }
    best.freq = 0
    for (i in 1:n.sims) {
      if(identical.groups && EqualPairs(answers[[i]]$Pairs, best.pairs$Pairs) || 
         !identical.groups && all(answers[[i]]$Pairs == best.pairs$Pairs)){
        best.freq = best.freq + 1
      }
    }
    answers <<- answers
    if(best.freq <= 1){
      if(best.freq == 1){
        print(paste("Only one simulation arrived at the best answer found, so there is a", 
                    "good chance there exists at least one other with a lower loss value"))
      }
      
      tryCatch({
        best.pairs$Pairs = AnalyzeAnswers(answers, loss.fn, sample.size, people, groups, n.groups, n.sims)
        best.pairs$Loss = loss.fn(best.pairs$Pairs, groups, n.groups, col.names)
      }, error = function(e){
        print("Unable to analyze answers to potentially find a better pairing")
      })
      
    } else if(best.freq > 1){
      print(paste(best.freq, 
                  "different simulations arrived at this answer, so there is a", 
                  "good chance it is the one that minimizes the loss function"))
    }
    print("Best answer:")
    print(best.pairs$Pairs)
    print(paste("Loss:", best.pairs$Loss))
    end.time = Sys.time()
    duration = difftime(end.time, start.time)
    print(paste(round(duration[[1]], 2), units(duration), "for", n.sims, "sims"))
    return(best.pairs$Pairs)
    
  }
  MakePairsAux = function(groups, people, loss.fn, prob.df, n.groups, sample.size, 
                          iter.level, identical.groups, alpha = 0.05, shuffled.ppl = NA, 
                          col.names = NA, progress.bar = TRUE){
    if(isTRUE(is.na(shuffled.ppl))){
      stop("No shuffled people")
    } else {
      if(iter.level == 1 && identical.groups){
        print(paste("First leader:", shuffled.ppl[1]))
      }
      global.shuffled.ppl <<- shuffled.ppl
    }
    # probs is the probability matrix of people being placed in a given group 
    # for this iteration, each entry is 1/groups
    print(paste0("Iteration ", iter.level, " (Sim ", sim.num, ")"))
    prob.df <<- prob.df
    
    if(isTRUE(is.na(col.names))){
      col.names = paste("Group", 1:n.groups)
    }
    nonzero.probs = CalcNonzeroProbs(prob.df, col.names)
    err = FALSE
    combo.upper.bound = Inf 
    if(err){
      if(identical.groups){
        print(("Restarting with new people shuffle"))
        new.sample.size = sample.size
      } else {
        print(("Restarting with doubled sample size"))
        new.sample.size = sample.size * 2
      }
      
      # default.prob.df
      return(MakePairsAux(groups, people, loss.fn, new.probs, n.groups, 
                          new.sample.size, 
                          1, identical.groups, alpha, shuffled.ppl = 
                            ifelse(rep(identical.groups, length(people)), 
                                   sample(people), people), col.names = col.names))
    }
    
    low.combos = combo.upper.bound < (n.groups * sample.size / 5)
    if(low.combos){
      iter.size = combo.upper.bound * 10
    } else {
      iter.size = round((n.groups + mean(nonzero.probs)) / 2 * sample.size)
    }
    print(paste("Size:", iter.size))
    
    
    
    tryCatch({
      sim.results = GeneratePairsData(groups, people, prob.df, n.groups, 
                                      iter.size, identical.groups, nonzero.probs, 
                                      progress.bar, shuffled.ppl = shuffled.ppl, 
                                      col.names = col.names)
      err = FALSE
    }, error = function(e){
      err <<- TRUE
      print(e)
    })
    if(err){
      print("Error generating pairs, restarting simulation with new shuffle")
      new.probs = default.prob.df
      return(MakePairsAux(groups, people, loss.fn, new.probs, n.groups, 
                          sample.size, 
                          1, identical.groups, alpha, shuffled.ppl = 
                            ifelse(rep(identical.groups, length(people)), 
                                   sample(people), people), col.names = col.names))
    }
    sim.data = OrganizeManyPairs(sim.results, people, col.names, iter.size, 
                                 loss.fn, iter.level, n.groups = n.groups)
    old.prob.df = prob.df
    distinct.pairs = sim.data %>% 
      mutate(Index = 1:nrow(sim.data)) %>% 
      distinct_at(vars(-Index), .keep_all = TRUE)
    if(nrow(distinct.pairs) < iter.size / 10){
      low.combos = TRUE
    } else if(identical.groups && nrow(distinct.pairs) < iter.size / 3){
      distinct.pairs = distinct.pairs %>% 
        arrange(Loss)
      distinct.loss = unique(distinct.pairs$Loss)
      if(length(distinct.loss) < iter.size / 10){
        failed = numeric()
        for (l in distinct.loss) {
          pairs.index = distinct.pairs[which(distinct.pairs$Loss == l), "Index"]
          passed = numeric()
          for (i in pairs.index) {
            new.pair = TRUE
            for (p in passed) {
              if(EqualPairs(sim.results[[i]], sim.results[[p]])){
                new.pair = FALSE
                break
              }
            }
            if(new.pair){
              passed = c(passed, i)
            } else {
              failed = c(failed, i)
            }
          }
        }
        distinct.pairs = distinct.pairs %>% 
          filter(!(Index %in% failed))
        if(nrow(distinct.pairs) < iter.size / 10){
          low.combos = TRUE
        }
      }
    }
    if(low.combos){
      print(paste("There were only", nrow(distinct.pairs), 
                  "permutations, so the one with the best loss was chosen"))
      new.probs = sim.results[[match(min(sim.data$Loss), sim.data$Loss)]]
    } else {

      new.probs = ProbabilityTests(people, n.groups, col.names, prob.df,
                                   sim.data, identical.groups, alpha)
      # overfull.groups = FindOverfullGroups(new.probs, groups, n.groups, col.names)
      # if(length(overfull.groups) > 0){
      #   new.probs = FixOverfullGroups(prob.df, new.probs, overfull.groups, sim.data,
      #                                 groups, n.groups, people, alpha, col.names)
      # }
    }
    tryCatch({
      if(all(unlist(new.probs == prob.df))){
        print("No progress made in iteration, generating more data")
        sim.data = GeneratePairsData(groups, people, prob.df, n.groups, 
                                     iter.size, identical.groups, nonzero.probs, 
                                     shuffled.ppl = shuffled.ppl, col.names = col.names) %>% 
          OrganizeManyPairs(people, col.names, iter.size, loss.fn, iter.level, 
                            progress.bar, n.groups) %>% 
          rbind(sim.data)
        new.probs = ProbabilityTests(people, n.groups, col.names, prob.df, sim.data, 
                                     identical.groups, alpha)
      }
    }, error = function(e){
      error.info <<- list("new.probs" = new.probs, "prob.df" = prob.df, 
                          "sim.data" = sim.data)
    })
    
    print(new.probs)
    set.people = people[which(new.probs == 1, arr.ind = TRUE)[,1]]
    print(paste("Set people:", paste(set.people, collapse = ", ")))
    new.set.ppl = set.people %>% 
      setdiff(people[which(old.prob.df == 1, arr.ind = TRUE)[,1]])
    print(paste("New set people:", paste(new.set.ppl, collapse = ", ")))
    elim.progress = sum(unlist(new.probs == 0)) / (length(people) * (n.groups - 1))
    print(paste("Progress:", length(set.people), "of", length(people), "people are set,",
                round(elim.progress, 3) * 100, "pct of groupings eliminated, best loss", signif(min(sim.data$Loss), 5)))
    if(any(unlist(prob.df == 0) & unlist(new.probs != 0))){
      debug.info = vector('list')
      debug.info$Sim.data = sim.data
      debug.info$Old.probs = prob.df
      debug.info$New.probs = new.probs
      print("Some zero entries became nonzero")
      return(debug.info)
    } else if(all(unlist(new.probs == prob.df))){
      if(identical.groups){
        print("No progress made in iteration, restarting with new people shuffle")
        new.sample.size = sample.size
      } else {
        print("No progress made in iteration, restarting with doubled sample size")
        new.sample.size = sample.size * 2
      }
      
      new.probs = default.prob.df
      return(MakePairsAux(groups, people, loss.fn, new.probs, n.groups, 
                          new.sample.size, 
                          1, identical.groups, alpha, shuffled.ppl = 
                            ifelse(rep(identical.groups, length(people)), 
                                   sample(people), people), col.names = col.names))
    } else if(all(unlist(new.probs[, col.names]) %in% c(0, 1))){
      return(new.probs)
    } else {
      # general case: recursive call
      return(MakePairsAux(groups, people, loss.fn, new.probs, n.groups, sample.size, 
                          iter.level + 1, identical.groups, alpha, 
                          shuffled.ppl = shuffled.ppl, col.names = col.names))
    }
  }
  AssignPairs = function(people, groups, n.groups, prob.df, nonzero.probs, attempt = 1, 
                         identical.groups, shuffled.ppl = NA, col.names){
    if(isTRUE(is.na(shuffled.ppl))){
      stop("No shuffled people")
    }
    # if(class(prob.df) == 'list'){
    #   return(AssignPairsBanList(people, groups, n.groups, prob.df, nonzero.probs))
    # }
    order.tiers = sort(unique(nonzero.probs))
    pick.order = character()
    for (i in order.tiers) {
      pick.order = c(pick.order, sample(people[which(nonzero.probs == i)]))
    }
    if(isTRUE(is.na(col.names))){
      col.names = paste("Group", 1:n.groups)
    }
    pairs.df = data.frame("Person" = people) %>% 
      AddCols(col.names, 0)
    capacities = groups
    tryCatch({
      for (person in pick.order) {
        person.index = match(person, people)
        probs = prob.df[person.index, col.names]
        ret.list = AssignOnePair(person, pairs.df, capacities, probs, col.names) 
        capacities = ret.list$Capacities
        pairs.df = ret.list$Pairs
      }
    }, error = function(e){
      if(attempt > 50){
        error.info <<- list("people" = people, 'groups' = groups, 'n.groups' = n.groups, 
                            'prob.df'= prob.df, 'nonzero.probs' = nonzero.probs, 
                            'identical.groups' = identical.groups, 
                            'pick.order' = pick.order, 'capacities' = capacities, 
                            'pairs.df' = pairs.df, 'message' = e$message)
        stop(paste("Pairs unsuccessful after 50 attempts\n", e$message))
      } else {
        if(attempt == 1){
          pick.errors <<- paste(person, match(person, pick.order))
        } else {
          pick.errors <<- c(pick.errors, paste(person, match(person, pick.order)))
        }
        pairs.df <<- AssignPairs(people, groups, n.groups, prob.df, nonzero.probs, attempt + 1, 
                                 identical.groups = identical.groups, 
                                 shuffled.ppl = shuffled.ppl, col.names = col.names)
      }
      
    })
    if(any(pairs.df != 0 & prob.df == 0)){
      print(pairs.df)
      print(prob.df)
      stop("Illegal pairing")
    }
    if(!LegalPairs(prob.df, pairs.df)){
      print(pairs.df)
      print(prob.df)
      stop("Illegal pairing")
    }
    
    # if(identical.groups){ 
    #   set.groups = which(prob.df == 1, arr.ind = TRUE)[,2] %>% unique() - 1
    #   free.groups = setdiff(1:n.groups, set.groups)
    #   group.sizes = numeric()
    #   lead.people = character()
    #   for (i in 1:n.groups) {
    #     group.index = match(people[which(pairs.df[, paste("Group", i)] == 1)], 
    #                         shuffled.ppl)
    #     group.sizes[i] = length(group.index)
    #     if(group.sizes[i] > 0){
    #       lead.people[i] = shuffled.ppl[min(group.index)]
    #     }
    #     
    #   }
    #   unq.group.sizes = unique(groups)
    #   new.groups = numeric()
    #   new.groups[set.groups] = set.groups
    #   # make it so if a person is locked into a group they stay there 
    #   for (j in unq.group.sizes) {
    #     # case where group is of size 0 and lead person is na
    #     group.nums = which(groups == j)
    #     free.new.groups = group.nums[which(group.nums %in% free.groups)]
    #     new.groups[free.new.groups] = free.new.groups[rank(match(lead.people[free.new.groups], 
    #                                                              shuffled.ppl))]
    #   }
    #   new.pairs = pairs.df
    #   for (i in 1:n.groups) {
    #     new.pairs[, paste("Group", new.groups[i])] = pairs.df[, paste("Group", i)]
    #   }
    #   prob.df <<- prob.df
    #   pairs.df <<- pairs.df
    #   if(LegalPairs(prob.df, new.pairs)){
    #     pairs.df = new.pairs
    #   } else {
    #     stop("Could not reorder groups")
    #   }
    # }
    return(pairs.df)
  }
  LegalPairs = function(prob.df, pairs.df){
    # make sure any 0 probs don't have a 1 pair
    return(!any(prob.df == 0 & pairs.df != 0))
  }
  AssignOnePair = function(person, pairs.df, capacities, probs, col.names){
    # person is a string representing the name of the person 
    # pairs.df is a data frame
    # capacities is an integer vector representing the groups' capacities
    # probs is a vector for the probability the person is placed in each group
    
    available.groups = which(capacities > 0)
    new.probs = probs[available.groups]
    if(sum(new.probs) == 0){
      stop(paste("All of", person, "potential groups are full"))
    }
    if(length(available.groups) > 1){
      new.probs = new.probs / sum(new.probs)
      assigned.group = sample(available.groups, 1, prob = new.probs)
    } else {
      assigned.group = available.groups
    }
    new.capacities = capacities
    new.capacities[assigned.group] = new.capacities[assigned.group] - 1
    
    person.index = match(person, pairs.df$Person)
    pairs.df[person.index, col.names[assigned.group]] = 1
    
    
    # return a list
    # first element is the updated pairs.df
    # second element is an updated (remaining) capacities vector
    # subtract one from the group this person was placed in 
    ret.list = vector('list', 2)
    ret.list$Pairs = pairs.df
    ret.list$Capacities = new.capacities
    return(ret.list)
  }
  GeneratePairsData = function(groups, people, prob.df, n.groups, iter.size, 
                               identical.groups, nonzero.probs, 
                               progress.bar = TRUE, shuffled.ppl = NA, col.names){
    if(isTRUE(is.na(shuffled.ppl))){
      stop("No shuffled people")
    }
    if(isTRUE(is.na(col.names))){
      col.names = paste("Group", 1:n.groups)
    }
    sim.results = vector('list', length = iter.size)
    
    if(progress.bar){
      print("Generating pairs data...")
      tryCatch({
        pb = txtProgressBar(min = 0, max = iter.size, style = 3)
      }, error = function(e){
        iter.size <<- iter.size
        isize.error.info <<- list('iter.size' = iter.size, 
                                  'nonzero.probs' = nonzero.probs)
        print(paste("Iter size", iter.size))
        stop(e)
      })
      
    }
    for (i in 1:iter.size) {
      sim.results[[i]] = AssignPairs(people, groups, n.groups, prob.df, 
                                     nonzero.probs, 1, identical.groups, 
                                     shuffled.ppl = shuffled.ppl, col.names = col.names)
      if(progress.bar){
        setTxtProgressBar(pb, i)
      }
    }
    if(progress.bar){
      close(pb)
    }
    return(sim.results)
  }
  OrganizePairs = function(sim.pairs, people, col.names, loss.fn, n.groups = NA){
    if(is.na(n.groups)){
      n.groups = length(col.names)
    }
    sim.row = data.frame()
    for (p in 1:length(people)) {
      sim.row[1, paste(people[p], "Group")] = 
        col.names[match(1, sim.pairs[p, col.names])]
    }
    sim.row$Loss = loss.fn(sim.pairs, groups, n.groups, col.names)
    return(sim.row)
  }
  OrganizeManyPairs = function(sim.results, people, col.names, iter.size, loss.fn, 
                               iter.level = NA, progress.bar = TRUE, n.groups = NA){
    if(is.na(n.groups)){
      n.groups = length(col.names)
    }
    if(iter.size != length(sim.results)){
      print("Iter.size changed to match length of sim.results")
      iter.size = length(sim.results)
    }
    if(progress.bar){
      print("Organizing pairs into a data frame and calculating loss...")
      pb = txtProgressBar(min = 0, max = iter.size, style = 3)
    }
    sim.data = data.frame()
    # first.row = OrganizePairs(sim.results[[i]], people, col.names, loss.fn, n.groups)
    # sim.data = matrix(nrow = iter.size, ncol = ncol(first.row)) %>% 
    #   as.data.frame() %>% 
    #     `colnames<-`(names(first.row))
    for (i in 1:iter.size) {
      sim.pairs = sim.results[[i]]
      sim.row = OrganizePairs(sim.pairs, people, col.names, loss.fn, n.groups)
      sim.data = sim.data %>% rbind(sim.row)
      # sim.data[i, ] = sim.row
      if(progress.bar){
        setTxtProgressBar(pb, i)
      }
    }
    if(progress.bar){
      close(pb)
    }
    sim.data$Iter = iter.level
    tryCatch({
      all.pairs.data <<- rbind(all.pairs.data, sim.data)
    }, error = function(e){
      all.pairs.data <<- sim.data
    })
    return(sim.data)
  }
  ProbabilityTests = function(people, n.groups, col.names, prob.df, sim.data, 
                              identical.groups, alpha = 0.05){
    new.probs = prob.df 
    group.levels = 1:n.groups
    if(identical.groups){
      set.groups = which(prob.df == 1, arr.ind = TRUE)[,2] %>% unique() - 1
    } else {
      set.groups = group.levels
    }
    free.groups = setdiff(group.levels, set.groups)
    
    for (p in 1:length(people)) {
      curr.probs = prob.df[p, col.names]
      if(all(curr.probs %in% c(0, 1))){
        new.probs[p, col.names] = curr.probs
      } else {
        group.col.name = paste(people[p], "Group")
        # mu0 = mean(sim data Loss)
        # xbar = mean(sim data Loss for specific group)
        # h0: mu = mu0 (where mu is derived from xbar)
        # ha: mu < mu0
        
        p.vals = numeric()
        for (g in set.groups) {
          p.vals[g] = SingleProbTest(sim.data$Loss, 
                                     which(sim.data[, group.col.name] == col.names[g]))
        }
        if(length(free.groups) > 0){
          p.vals[(g+1):n.groups] = length(free.groups) * 
            SingleProbTest(sim.data$Loss,
                           which(sim.data[, group.col.name] %in% free.groups))
        } # test all free groups as one 
        new.probs[p, col.names] = ifelse(p.vals < alpha, 0, 1) %>% 
          ProbVector()
      }
    }
    return(new.probs)
  }
  SingleProbTest = function(loss.data, test.index){
    return(ifelse(length(test.index) == 0, 0, 
                  tryCatch({
                    t.test(loss.data[test.index], loss.data[-test.index], 
                           alternative = 'greater')$p.value
                  }, error = function(e){
                    return(0.5)
                  }))
    )
  }
  CalcNonzeroProbs = function(prob.df, col.names = NA){
    if(isTRUE(is.na(col.names))){
      col.names = names(prob.df[-1])
    }
    nonzero.probs = numeric()
    for (i in 1:nrow(prob.df)) {
      nonzero.probs[i] = sum(prob.df[i, col.names] != 0)
    }
    return(nonzero.probs)
  }
  FindOverfullGroups = function(prob.df, groups, n.groups, col.names){
    overfull.groups = numeric()
    for (i in 1:n.groups) {
      if(sum(prob.df[, col.names[i]] == 1) > groups[i]){
        overfull.groups = c(overfull.groups, i)
      }
    }
    return(overfull.groups)
  }
} # functions
# preferences file 
folder.path = "/Users/malexk999/Desktop/Cloud desktop/University of Maryland/Frat/Big little 2021/Algorithm explanation"
file.name = "Example data.xlsx"


twins.elig = c()
# twins.elig = "all"
# twins.elig = c("Howard", "George")

# make it so that littles MUST be paired with a big they requested
use.default.groups.shortcut = FALSE # TRUE or FALSE

{
  start.time = Sys.time()
  library(readxl)
  library(dplyr)
  FixUnrankedBrother = function(choice.index){
    return(ifelse(is.na(choice.index), brotherUnranked, choice.index))
  }
  FixUnrankedPledge = function(choice.index){
    return(ifelse(is.na(choice.index), pledgeUnranked, choice.index))
  }
  # unrankedValue = 5
  sheet.path = paste0(folder.path, "/", file.name)
  brother.data = as.data.frame(read_xlsx(sheet.path, sheet = "Brother"))
  pledge.data = as.data.frame(read_xlsx(sheet.path, sheet = "Pledge"))
  nPledges = nrow(pledge.data)
  pledgeChoices = ncol(pledge.data) - 1 # pledges choosing brothers
  nBrothers = nrow(brother.data)
  brotherChoices = ncol(brother.data) - 1 # brothers choosing pledges
  # value of a brother getting a pledge he didn't choose 
  pledgeUnranked = (brotherChoices + nPledges) / 2
  brotherUnranked = (pledgeChoices + nBrothers) / 2
  people = pledge.data$Pledge # lol pledges aren't people
  brother.names = brother.data$Brother
  if(isTRUE(twins.elig == "all")){
    groups = rep(2, nBrothers)
  } else {
    groups = rep(1, nBrothers)
    groups[match(twins.elig, brother.names)] = 2
  } # who can get twins? 
  if(use.default.groups.shortcut){
    default.groups = vector('list')
    for (i in 1:nrow(pledge.data)) {
      default.groups[[people[i]]] = match(pledge.data[i, -1], brother.names)
    }
  } else {
    default.groups = NULL
  }
  
  loss.fn = function(pairs.df, groups, n.groups, col.names = names(pairs.df)[-1]){
    total.loss = 0
    for (i in 1:n.groups) {
      group.members = pairs.df %>% 
        filter_at(col.names[i], function(x) x > 0) %>% 
        select(Person) %>% 
        unlist() %>% 
        as.character()
      brother.requests = brother.data[i, -1] 
      for (person in group.members) {
        pledge.requests = pledge.data[match(person, people), -1] %>% as.character() 
        big.choice = match(brother.names[i], pledge.requests) %>% 
          FixUnrankedBrother()
        little.choice = match(person, brother.requests) %>% FixUnrankedPledge()
        total.loss = total.loss + big.choice ^ 2 + little.choice ^ 2
      }
    }
    return(total.loss)
  }
  {
    called.fns = NULL
    all.pairs.data = data.frame()
    answer = MakePairs(people, groups, loss.fn, default.groups = default.groups, 
                       group.names = brother.names)
    first.iter.data = all.pairs.data %>% 
      filter(Iter == 1) %>% 
      select(-Iter) %>% 
      distinct()
    all.pairs.data = all.pairs.data %>% 
      select(-Iter) %>%
      distinct()
    # print(called.fns)
    end.time = Sys.time()
    time.diff = end.time - start.time
    print(time.diff)
    beepr::beep()
  } # where the magic happens 
  {
    # unranked value = (5 + 3) / 2 = 4 for brothers and pledges
    
    # Samuel / Garrett (Samuel picked Garrett 3rd, Garrett had Samuel 2nd) (3, 2)
    # Jacob / Chad (3, 1)
    # Howard / Andrew (4, 2) (Howard did not choose Andrew: value is 4)
    # George / Tyler  (3, 2)
    # Richard / Gabe (1, 3)
    # 3^2 + 2^2 + 3^2 + 1^2 + 4^2 + 2^2 + 3^2 + 2^2 + 1^2 + 3^2 = 66
    example.pairing = data.frame("Person" = people) %>% 
      AddCols(paste("Group", 1:5), 0) %>% 
      mutate("Group 1" = as.numeric(row_number() == 3)) %>% 
      mutate("Group 2" = as.numeric(row_number() == 2)) %>% 
      mutate("Group 3" = as.numeric(row_number() == 1)) %>% 
      mutate("Group 4" = as.numeric(row_number() == 4)) %>% 
      mutate("Group 5" = as.numeric(row_number() == 5))
    # loss.fn(example.pairing, groups, ncol(example.pairing) - 1, brother.names) # 66
    
    
  } # example 
} # where the magic happens 


{
  squared.loss = function(pairs.df, groups, n.groups, col.names = names(pairs.df)[-1]){
    total.loss = 0
    for (i in 1:n.groups) {
      group.members = pairs.df %>% 
        filter_at(col.names[i], function(x) x > 0) %>% 
        select(Person) %>% 
        unlist() %>% 
        as.character()
      brother.requests = brother.data[i, -1] 
      for (person in group.members) {
        pledge.requests = pledge.data[match(person, people), -1] %>% as.character() 
        big.choice = match(brother.names[i], pledge.requests) %>% 
          FixUnrankedBrother()
        little.choice = match(person, brother.requests) %>% FixUnrankedPledge()
        total.loss = total.loss + big.choice ^ 2 + little.choice ^ 2
      }
    }
    return(total.loss)
  }
  linear.loss = function(pairs.df, groups, n.groups, col.names = names(pairs.df)[-1]){
    total.loss = 0
    for (i in 1:n.groups) {
      group.members = pairs.df %>% 
        filter_at(col.names[i], function(x) x > 0) %>% 
        select(Person) %>% 
        unlist() %>% 
        as.character()
      brother.requests = brother.data[i, -1] 
      for (person in group.members) {
        pledge.requests = pledge.data[match(person, people), -1] %>% as.character() 
        big.choice = match(brother.names[i], pledge.requests) %>% 
          FixUnrankedBrother()
        little.choice = match(person, brother.requests) %>% FixUnrankedPledge()
        total.loss = total.loss + big.choice ^ 1 + little.choice ^ 1
      }
    }
    return(total.loss)
  }
  sqrt.loss = function(pairs.df, groups, n.groups, col.names = names(pairs.df)[-1]){
    total.loss = 0
    for (i in 1:n.groups) {
      group.members = pairs.df %>% 
        filter_at(col.names[i], function(x) x > 0) %>% 
        select(Person) %>% 
        unlist() %>% 
        as.character()
      brother.requests = brother.data[i, -1] 
      for (person in group.members) {
        pledge.requests = pledge.data[match(person, people), -1] %>% as.character() 
        big.choice = match(brother.names[i], pledge.requests) %>% 
          FixUnrankedBrother()
        little.choice = match(person, brother.requests) %>% FixUnrankedPledge()
        total.loss = total.loss + big.choice ^ 0.5 + little.choice ^ 0.5
      }
    }
    return(total.loss)
  }
  ManualPair = function(file.name){
    require(xlsx)
    require(readxl)
    folder.path = "/Users/malexk999/Desktop/Cloud desktop/University of Maryland/Frat/Big little 2021/Possible pairings"
    read.path = paste0(folder.path, "/Manual in/", file.name, ".xlsx")
    manual.in = as.data.frame(read_excel(read.path))
    pairs.df = default.prob.df
    for(p in people){
      i = match(p, manual.in$Little)
      g = match(manual.in[i, "Big"], brother.names)
      pairs.df[i, -1] = 0
      pairs.df[i, g + 1] = 1
    } # build pairs df 
    pledge.picks = numeric(length = length(people))
    brother.picks = numeric(length = length(people))
    bigs = numeric(length = length(brother.names))
    for(i in 1:nrow(pairs.df)){
      bigs[i] = brother.names[match(1, answer[i, -1])]
      pledge.picks[i] = match(bigs[i], pledge.data[i, 2:6] %>% as.character())
      brother.picks[i] = match(people[i], brother.data[match(
        bigs[i], brother.names), 2:6] %>% as.character())
    }
    manual.in$Big.pick = brother.picks
    manual.in$Little.pick = pledge.picks
    loss.types = c("Squared", "Linear", "Sqrt")
    loss.vals = c(squared.loss(pairs.df, groups, length(groups)), 
                  linear.loss(pairs.df, groups, length(groups)), 
                  sqrt.loss(pairs.df, groups, length(groups)))
    loss.sheet = data.frame("Type" = loss.types, 
                            "Loss" = loss.vals)
    out.path = paste0(folder.path, "/Manual out/", file.name, ".out.xlsx")
    write.xlsx(manual.in, out.path, sheetName = "Pairings", row.names = FALSE)
    suppressWarnings(write.xlsx(loss.sheet, out.path, sheetName = "Loss", append = TRUE, row.names = FALSE))
    print(loss.sheet)
  }
  # best found solution loss: 314, 72, 40.54534
} # advanced


{
  # folder.path = "/Users/malexk999/Desktop/Cloud desktop/University of Maryland/Frat/Big little 2021"
  # file.name = "big little 2021.xlsx"
  # twins.elig = "Taylor Bennett"
  
  # use default groups thing 
  # 1 minute for TRUE, 2 minutes for FALSE
  # this did not work for big little spring 2021 data 
  
} # old comments


