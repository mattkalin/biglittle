# This algorithm was created by Matt Kalin in April 2020 
# It was used to create big/little pairings for Sigma Alpha Mu fraternity at 
# the University of Maryland, College Park during the Spring 2021 semester 
# Contact: kalinmatt4@gmail.com 
# git repository: https://github.com/mattkalin/biglittle


{
  called.fns = NULL
  ProbVector = function(x){
    called.fns <<- union(called.fns, "ProbVector")
    return(x / sum(x))
  }
  AddCols = function(df, new.cols, default.val = NA){
    called.fns <<- union(called.fns, "AddCols")
    df[, new.cols] = default.val
    return(df)
  }
  MakePairs = function(people, groups, loss.fn, sample.size = 50,
                       n.sims = 5, identical.groups = FALSE, 
                       default.groups = NULL, alpha = 0.05){
    called.fns <<- union(called.fns, "MakePairs")
    
    # people is a character vector of everyone's names 
    # loss.fn is the loss function (high output is bad, low output is good)
    # capacities is EITHER 
    # a vector of the same length as groups that represents the capacity of each group OR
    # an integer representing the number of groups
    # first.first: if TRUE, A will always be placed in group 1
    
    # consider making iter.size = 50 (or other constant) * n.groups 
    start.time <<- Sys.time()
    library(dplyr)
    n.groups = length(groups)
    # iter.size = sample.size * n.groups
    if(n.groups == 1){
      n.groups = groups
      groups = rep(ceiling(length(people)/n.groups), n.groups)
    } # else {
    #   n.groups = length(groups)
    # }
    # if(identical.groups & is.null(default.groups)){
    #   prob.df = vector('list', length = length(people))
    #   names(prob.df) = people
    #   # ban list 
    # } else {
    col.names = paste("Group", 1:n.groups)
    prob.df = data.frame("Person" = people) %>% 
      AddCols(col.names, 1/n.groups)
    # }
    if(!is.null(default.groups)){
      for (p in names(default.groups)) {
        p.grp = default.groups[[p]]
        p.index = match(p, people)
        prob.df[p.index, col.names] = 0
        prob.df[p.index, paste("Group", p.grp)] = 1 / length(p.grp)
      }
      # for (i in 1:nrow(default.groups)) {
      #   row.index = match(default.groups[i, "Person"], people)
      #   prob.df[row.index, col.names] = 0
      #   prob.df[row.index, paste("Group", default.groups[i, "Group"])] = 1
      # }
    }
    default.prob.df <<- prob.df
    # consider making it weighted based on the groups' capacities 
    # ex: more likely to be placed in bigger group 
    # although I think that's gonna happen anyway tbh 
    
    # if(first.first){
    #   prob.df[1, col.names] = c(1, rep(0, n.groups - 1))
    # }
    # top.total = iter.size * top.prop # the total number of sims in the top proportion
    answers = vector('list', n.sims)
    best.pairs = vector('list', 2)
    best.pairs$Loss = Inf
    for (k in 1:n.sims) {
      sim.num <<- k
      # answers[[k]]$Pairs = guess.answer
      print(paste("Simulation", k))
      answers[[k]]$Pairs = # tryCatch({
        MakePairsAux(groups, people, loss.fn, prob.df, n.groups, sample.size, 
                     # top.total, conv.diff, 
                     1, 
                     # first.first, method, 
                     identical.groups, alpha, shuffled.ppl = 
                       ifelse(rep(identical.groups, length(people)), 
                              sample(people), people))
      # }, error = function(e){
      # print("Error in MakePairsAux, trying again")
      # return(MakePairsAux(groups, people, loss.fn, prob.df, n.groups, sample.size, 
      #                     # top.total, conv.diff, 
      #                     1, 
      #                     # first.first, method, 
      #                     identical.groups))
      #   stop(e)
      # })
      
      if(!is.null(answers[[k]]$Pairs$Sim.data)){
        return(answers[[k]]$Pairs)
      }
      answers[[k]]$Loss = loss.fn(answers[[k]]$Pairs, groups, n.groups) 
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
        new.pairs[i, paste("Group", best.row[, i])] = 1
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
    # answer = best.pairs$Pairs
    if(best.freq <= 1){
      if(best.freq == 1){
        print(paste("Only one simulation arrived at the best answer found, so there is a", 
                    "good chance there exists at least one other with a lower loss value"))
      }
      
      tryCatch({
        best.pairs$Pairs = AnalyzeAnswers(answers, loss.fn, sample.size, people, groups, n.groups, n.sims)
        best.pairs$Loss = loss.fn(best.pairs$Pairs, groups, n.groups)
      }, error = function(e){
        print("Unable to analyze answers to potentially find a better pairing")
      })
      
    } else if(best.freq > 1){
      print(paste(best.freq, 
                  "different simulations arrived at this answer, so there is a", 
                  "good chance it is the one that minimizes the loss function"))
      # return(best.pairs$Pairs)
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
                          iter.level, identical.groups, alpha = 0.05, shuffled.ppl = NA){
    called.fns <<- union(called.fns, "MakePairsAux")
    if(isTRUE(is.na(shuffled.ppl))){
      stop("No shuffled people")
    } else {
      if(iter.level == 1 && identical.groups){
        print(paste("First leader:", shuffled.ppl[1]))
      }
      global.shuffled.ppl <<- shuffled.ppl
    }
    # probs is the probability matrix of people being placed in a given group 
    # for the first iteration, each entry is 1/groups
    print(paste0("Iteration ", iter.level, " (Sim ", sim.num, ")"))
    prob.df <<- prob.df
    
    # preserve: nonzero.probs (argument), loss.results (do later)
    # iter.size (replace with sample.size in argumetns)
    
    # START REPLACEMENT FUNCTION HERE
    # args: groups, people, prob.df, n.groups, sample.size, identical.groups
    # do loss function calculations when organizing the sim.data df 
    col.names = paste("Group", 1:n.groups)
    # pairs.totals = data.frame("Person" = people) %>% 
    #   AddCols(col.names, 0)
    # sim.results = vector('list', length = iter.size)
    # loss.results = vector('numeric', length = iter.size)
    if(class(prob.df) == 'list'){
      nonzero.probs = numeric()
      # prob.df is ban.list, nonzero.probs is team.bans
      for (i in 1:length(prob.df)) {
        nonzero.probs[i] = length(prob.df[[i]])
      }
      iter.size = round(sample.size * 
                          (length(people) - 1 - mean(nonzero.probs) / 2) / mean(groups - 1))
    } else {
      nonzero.probs = CalcNonzeroProbs(prob.df, col.names)
      err = FALSE
      # tryCatch({
      #   combo.upper.bound = EstimateCombinations(prob.df, groups)
      # }, error = function(e){
      #   err <<- TRUE
      #   print(e$message)
      # })
      combo.upper.bound = Inf # trying a different approach 
      if(err){
        if(identical.groups){
          print(("Restarting with new people shuffle"))
          new.sample.size = sample.size
        } else {
          print(("Restarting with doubled sample size"))
          new.sample.size = sample.size * 2
        }
        
        # maybe just continue where it is now instead of completely restarting 
        # new.probs = data.frame("Person" = people) %>% 
        #   AddCols(col.names, 1/n.groups)
        default.prob.df
        return(MakePairsAux(groups, people, loss.fn, new.probs, n.groups, 
                            new.sample.size, 
                            # top.total, conv.diff, 
                            1, 
                            # first.first, 2, 
                            identical.groups, alpha, shuffled.ppl = 
                              ifelse(rep(identical.groups, length(people)), 
                                     sample(people), people)))
      }
      
      low.combos = combo.upper.bound < (n.groups * sample.size / 5)
      if(low.combos){
        iter.size = combo.upper.bound * 10
      } else {
        iter.size = round((n.groups + mean(nonzero.probs)) / 2 * sample.size)
      }
    } # nonzero probabilities
    # iter.size = round((2 * n.groups - sum(prob.df == 0)/length(people)) / 2) * sample.size
    print(paste("Size:", iter.size))
    
    
    # for (i in 1:iter.size) {
    #   sim.pairs = AssignPairs(people, groups, n.groups, prob.df, nonzero.probs, 
    #                           identical.groups = identical.groups)
    #   loss.results[i] = loss.fn(sim.pairs, groups, n.groups)
    #   sim.results[[i]] = sim.pairs
    # }
    tryCatch({
      sim.results = GeneratePairsData(groups, people, prob.df, n.groups, 
                                      iter.size, identical.groups, nonzero.probs, 
                                      shuffled.ppl = shuffled.ppl)
      err = FALSE
    }, error = function(e){
      err <<- TRUE
      print(e)
    })
    if(err){
      print("Error generating pairs, restarting simulation with new shuffle")
      # new.probs = data.frame("Person" = people) %>% 
      #   AddCols(col.names, 1/n.groups)
      new.probs = default.prob.df
      return(MakePairsAux(groups, people, loss.fn, new.probs, n.groups, 
                          sample.size, 
                          # top.total, conv.diff, 
                          1, 
                          # first.first, 2, 
                          identical.groups, alpha, shuffled.ppl = 
                            ifelse(rep(identical.groups, length(people)), 
                                   sample(people), people)))
    }
    # if(method == 1){
    # recursive simulation (original method)
    
    # loss.ranks = rank(loss.results)
    # top.pairs = sim.results[which(loss.ranks <= top.total)]
    # if(length(top.pairs) == 0){
    #   top.index = which(loss.ranks <= top.total)
    #   for (i in 1:iter.size) {
    #     if(loss.ranks[i] <= top.total){
    #       top.pairs = c(top.pairs, sim.results[[i]])
    #     }
    #   }
    # }
    # # if(length(top.pairs) < top.total / 2){
    # #   loss.results <<- loss.results
    # #   sim.results <<- sim.results
    # #   top.pairs <<- top.pairs
    # #   print(paste("Top pairs length",length(top.pairs)))
    # #   print(loss.ranks)
    # #   stop(paste("top.pairs of length", length(top.pairs)))
    # #   # return(MakePairsAux(groups, people, loss.fn, prob.df, n.groups, iter.size, 
    # #   #                     top.prop, conv.diff, iter.level))
    # # } else {
    # top.total = length(top.pairs)
    # # }
    # for (i in 1:top.total) {
    #   # tryCatch({
    #   pairs.totals[, col.names] = pairs.totals[, col.names] + top.pairs[[i]][, col.names]
    #   # }, error = function(e){
    #   #   print(paste("i =", i))
    #   #   print(paste("top.total =", top.total))
    #   #   error <<- e
    #   #   print(e)
    #   # })
    #   
    # }
    # new.probs = pairs.totals %>% 
    #   mutate_at(col.names, function(x) x / top.total) # this should work but not 100% sure 
    # prob.diff = new.probs[, col.names] - prob.df[, col.names]
    # if(any(abs(prob.diff) > conv.diff)){
    #   # abs.prob.diff = abs(unlist(prob.diff))
    #   # print(paste("Mean diff =", mean(abs.prob.diff)))
    #   # print(paste("Max diff =", max(abs.prob.diff)))
    #   all.prob.diff[[iter.level]] <<- prob.diff
    #   return(MakePairsAux(groups, people, loss.fn, new.probs, n.groups, sample.size, 
    #                       top.total, conv.diff, iter.level + 1, first.first))
    # } else {
    #   # also have it stop if each row has a 1 and 0s 
    #   return(new.probs)
    # }
    # } else if(method == 2){
    sim.data = OrganizeManyPairs(sim.results, people, col.names, iter.size, 
                                 loss.fn, iter.level, n.groups = n.groups)
    # sim.data = data.frame()
    # for (i in 1:iter.size) {
    #   sim.pairs = sim.results[[i]]
    #   sim.row = data.frame()
    #   for (p in 1:length(people)) {
    #     sim.row[1, paste(people[p], "Group")] = match(1, sim.pairs[p, col.names])
    #   }
    #   sim.row$Loss = loss.fn(sim.pairs, groups, n.groups)
    #   sim.data = sim.data %>% rbind(sim.row)
    # }
    # sim.data$Iter = iter.level
    # tryCatch({
    #   all.pairs.data <<- rbind(all.pairs.data, sim.data)
    # }, error = function(e){
    #   all.pairs.data <<- sim.data
    # })
    
    # model.people = head(people, -1)
    # if(first.first){
    #   model.people = model.people[-1]
    # } 
    # model.formula = paste("Loss ~", paste0("factor(`", model.people, " Group`)",
    #                                        collapse = " + ")) %>% 
    #   as.formula()
    # m = lm(model.formula, data = sim.data)
    # linear regression approach didn't work 
    
    # group.levels = as.character(1:n.groups)
    # for (p in 1:length(people)) {
    #   curr.probs = prob.df[p, col.names]
    #   if(all(curr.probs %in% c(0, 1))){
    #     new.probs[p, col.names] = curr.probs
    #   } else {
    #     group.col.name = paste(people[p], "Group")
    #     sim.data[, group.col.name] = sim.data[, group.col.name] %>% 
    #       as.factor()
    #     # group.levels = levels(sim.data[, group.col.name])
    #     loss.given.group = vector('list', length = length(group.levels))
    #     # unpooled.variance = 0
    #     # sample.means = numeric()
    #     
    #     # mu0 = mean(sim data Loss)
    #     # xbar = mean(sim data Loss for specific group)
    #     # h0: mu = mu0 (where mu is derived from xbar)
    #     # ha: mu < mu0
    #     
    #     p.vals = numeric()
    #     for (g in group.levels) {
    #       loss.given.group[[g]] = sim.data %>% 
    #         filter_at(group.col.name, function(x) x == g) %>% 
    #         select(Loss) %>% 
    #         unlist()
    #       # unpooled.variance = unpooled.variance + 
    #       #   sd(loss.given.group[[g]]) ^ 2 / length(loss.given.group[[g]])
    #       # sample.means[g] = mean(loss.given.group[[g]])
    #       
    #       # eilminate pairings that associate with worse than average loss values 
    #       if(length(loss.given.group[[g]]) < 1){
    #         p.vals[[g]] = 0
    #       } else {
    #         tryCatch({
    #           test.result = t.test(loss.given.group[[g]], mu = mean(sim.data$Loss), 
    #                                alternative = 'greater')
    #           p.vals[[g]] = test.result$p.value
    #         }, error = function(e){
    #           p.vals[[g]] <<- 0.5
    #         })
    #         
    #       }
    #       # if(p.vals[[g]] > pValueThreshold){
    #       #   p.vals[[g]] = 1
    #       # }
    #     }
    #     # new.probs[p, col.names] = (1 - p.vals) / sum(1 - p.vals)
    #     new.probs[p, col.names] = ifelse(p.vals < (0.05), 0, 1) %>% 
    #       ProbVector()
    #   }
    # }
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
      if(identical.groups){
        if(class(prob.df) == 'data.frame'){
          new.data = tryCatch({
            SetGroups(people, n.groups, col.names, prob.df, sim.data, 
                      alpha, shuffled.ppl = shuffled.ppl)
          }, error = function(e){
            error.info <<- list('prob.df' = prob.df, 'sim.data' = sim.data, 
                                'shuffled.ppl' = shuffled.ppl)
          })
          
          prob.df = new.data$Prob.df
          sim.data = new.data$Sim.data
        } 
        
        # remove identical groups stuff from other parts of algorithm 
      }
      if(class(prob.df) == 'data.frame'){
        new.probs = ProbabilityTests(people, n.groups, col.names, prob.df, 
                                     sim.data, identical.groups, alpha)
        overfull.groups = FindOverfullGroups(new.probs, groups, n.groups)
        if(length(overfull.groups) > 0){
          # x <<- list("New.Probs" = new.probs, 
          #         "Old.Probs" = prob.df, 
          #         "Sim.Data" = sim.data)
          # stop("Info saved in x")
          new.probs = FixOverfullGroups(prob.df, new.probs, overfull.groups, sim.data, 
                                        groups, n.groups, people, alpha)
        }
      } else {
        new.ban.list = PeoplePairTests(sim.data, people, prob.df, alpha)
        global.ban.list <<- new.ban.list
      }
      
    }
    if(class(prob.df) == 'list'){
      ban.list = prob.df
      if(ListsEqual(ban.list, new.ban.list)){
        print("No progress made in iteration, generating more data")
        sim.data = GeneratePairsData(groups, people, prob.df, n.groups, 
                                     iter.size, identical.groups, nonzero.probs) %>% 
          OrganizeManyPairs(people, col.names, iter.size, loss.fn, iter.level) %>% 
          rbind(sim.data)
        new.ban.list = PeoplePairTests(people, n.groups, col.names, prob.df, sim.data, alpha)
      }
      # PrintList(new.ban.list)
      # I should prob print something tbh, but idk
      if(ListsEqual(ban.list, new.ban.list)){
        print("No progress made in iteration, restarting with doubled sample size")
        return(MakePairsAux(groups, people, loss.fn, list(), n.groups, 
                            sample.size * 2, 
                            # top.total, conv.diff, 
                            1, 
                            # first.first, 2, 
                            identical.groups, alpha))
      }
      group.leaders = FindCliques(new.ban.list, n.groups, people, print.msg = TRUE)
      
      if(length(group.leaders) == 0){
        return(MakePairsAux(groups, people, loss.fn, new.ban.list, n.groups, sample.size, 
                            # top.total, conv.diff, 
                            iter.level + 1, 
                            # first.first, 2, 
                            identical.groups, alpha))
      } else {
        group.leaders = group.leaders[[1]] # don't think it matters much 
        prob.df = BuildProbDfFromBanList(ban.list, group.leaders, 
                                         n.groups, people, col.names)
        return(MakePairsAux(groups, people, loss.fn, prob.df, n.groups, sample.size, 
                            # top.total, conv.diff, 
                            iter.level + 1, 
                            # first.first, 2, 
                            identical.groups, alpha))
      }
      # finish this 
    }
    tryCatch({
      if(all(unlist(new.probs == prob.df))){
        print("No progress made in iteration, generating more data")
        sim.data = GeneratePairsData(groups, people, prob.df, n.groups, 
                                     iter.size, identical.groups, nonzero.probs, 
                                     shuffled.ppl = shuffled.ppl) %>% 
          OrganizeManyPairs(people, col.names, iter.size, loss.fn, iter.level) %>% 
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
      
      # maybe just continue where it is now instead of completely restarting 
      # new.probs = data.frame("Person" = people) %>% 
      #   AddCols(col.names, 1/n.groups)
      new.probs = default.prob.df
      return(MakePairsAux(groups, people, loss.fn, new.probs, n.groups, 
                          new.sample.size, 
                          # top.total, conv.diff, 
                          1, 
                          # first.first, 2, 
                          identical.groups, alpha, shuffled.ppl = 
                            ifelse(rep(identical.groups, length(people)), 
                                   sample(people), people)))
      # return(MakePairs(people, groups, loss.fn, sample.size * 2, 
      #                  conv.diff = conv.diff, first.first = first.first, method = 2, 
      #                  n.sims = 1))
    } else if(all(unlist(new.probs[, col.names]) %in% c(0, 1))){
      return(new.probs)
    } else {
      return(MakePairsAux(groups, people, loss.fn, new.probs, n.groups, sample.size, 
                          # top.total, conv.diff, 
                          iter.level + 1, 
                          # first.first, 2, 
                          identical.groups, alpha, shuffled.ppl = shuffled.ppl))
    }
    # }
    
  }
  AssignPairs = function(people, groups, n.groups, prob.df, nonzero.probs, attempt = 1, 
                         identical.groups, shuffled.ppl = NA){
    called.fns <<- union(called.fns, "AssignPairs")
    # pick.order = sample(people, length(people))
    # if(first.first){
    #   # check if any rows in prob.df are all 1s and 0s and put those first 
    #   # pass it in as an argument  
    #   pick.order = pick.order[pick.order != people[1]] %>% 
    #     append(people[1], after = 0)
    # }
    if(isTRUE(is.na(shuffled.ppl))){
      stop("No shuffled people")
    }
    if(class(prob.df) == 'list'){
      return(AssignPairsBanList(people, groups, n.groups, prob.df, nonzero.probs))
    }
    order.tiers = sort(unique(nonzero.probs))
    pick.order = character()
    for (i in order.tiers) {
      pick.order = c(pick.order, sample(people[which(nonzero.probs == i)]))
    }
    col.names = paste("Group", 1:n.groups)
    pairs.df = data.frame("Person" = people) %>% 
      AddCols(col.names, 0)
    capacities = groups
    tryCatch({
      for (person in pick.order) {
        person.index = match(person, people)
        probs = prob.df[person.index, col.names]
        ret.list = AssignOnePair(person, pairs.df, capacities, probs) 
        capacities = ret.list$Capacities
        pairs.df = ret.list$Pairs
        # if(pick.order[1] == "Jonas"){
        #   stop("Jonas is first, this is just a debugging tool")
        # }
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
                                 shuffled.ppl = shuffled.ppl)
      }
      
    })
    if(any(pairs.df != 0 & prob.df == 0)){
      print(pairs.df)
      print(prob.df)
      stop("Illegal pairing")
    }
    # if(any(pairs.df[10, 2:3] != 0)){
    #   print(pairs.df)
    #   print(prob.df)
    #   print("Illegal pairing for DEN")
    #   stop("Illegal pairing for DEN")
    # }
    if(!LegalPairs(prob.df, pairs.df)){
      print(pairs.df)
      print(prob.df)
      stop("Illegal pairing")
    }
    # bug: 
    # example: H can go to group 3 but not group 2
    # H is placed in group 3
    # group 2 and 3 are swtiched in the identical.groups step
    # H is now in group 2
    # potential fixes: 
    # 1
    # do itentical.groups in AssignOnePair, and if group 2 is empty, put there instead of 3
    # drawback: need to rework the pick.order step 
    # 2
    
    if(identical.groups){ # PREVIOUSLY: if(identical.groups)
      # THIS ONLY APPLIES TO GROUPS OF SAME SIZE 
      # lead person is the person in each group listed first in 'people' 
      
      set.groups = which(prob.df == 1, arr.ind = TRUE)[,2] %>% unique() - 1
      free.groups = setdiff(1:n.groups, set.groups)
      group.sizes = numeric()
      lead.people = character()
      for (i in 1:n.groups) {
        group.index = match(people[which(pairs.df[, paste("Group", i)] == 1)], 
                            shuffled.ppl)
        group.sizes[i] = length(group.index)
        if(group.sizes[i] > 0){
          lead.people[i] = shuffled.ppl[min(group.index)]
        }
        
      }
      unq.group.sizes = unique(groups)
      new.groups = numeric()
      new.groups[set.groups] = set.groups
      # make it so if a team is locked into a group they stay there 
      for (j in unq.group.sizes) {
        # case where group is of size 0 and lead person is na
        group.nums = which(groups == j)
        free.new.groups = group.nums[which(group.nums %in% free.groups)]
        new.groups[free.new.groups] = free.new.groups[rank(match(lead.people[free.new.groups], 
                                                                 shuffled.ppl))]
      }
      new.pairs = pairs.df
      for (i in 1:n.groups) {
        new.pairs[, paste("Group", new.groups[i])] = pairs.df[, paste("Group", i)]
      }
      prob.df <<- prob.df
      pairs.df <<- pairs.df
      if(LegalPairs(prob.df, new.pairs)){
        pairs.df = new.pairs
      } else {
        stop("Could not reorder groups")
        # adj.matrix = matrix(nrow = n.groups, ncol = n.groups, data = 0)
        # eligible.groups = list()
        # for (i in 1:n.groups) {
        #   member.rows = prob.df[which(pairs.df[, paste("Group", i)] == 1), ]
        #   eligible.groups[[i]] = numeric()
        #   for (j in 1:n.groups) {
        #     if(j != i && all(member.rows[, paste("Group", j)] > 0)){
        #       adj.matrix[i, j] = 1
        #       eligible.groups[[i]] = c(eligible.groups[[i]], j)
        #       # i is src, j is dest 
        #     }
        #   }
        #   # eligible.groups[[i]] = sort(eligible.groups[[i]], method = 'bubble')
        # }
        # 
      }
    }
    return(pairs.df)
  }
  LegalPairs = function(prob.df, pairs.df){
    called.fns <<- union(called.fns, "LegalPairs")
    # make sure any 0 probs don't have a 1 pair
    return(!any(prob.df == 0 & pairs.df != 0))
  }
  AssignOnePair = function(person, pairs.df, capacities, probs){
    called.fns <<- union(called.fns, "AssignOnePair")
    # person is a string representing the name of the person [REMOVED]
    # person.index is the person's index in the person vector 
    # number is the group number the person is being assigned to [REMOVED]
    # pairs.df is a data frame
    # first column is a list of all the people's names
    # next columns are binary: 1 if person is in group, 0 if they're in another group
    # NA if not yet assigned 
    # capacities is an integer vector representing the groups' capacities
    # probs is a vector for the probability the person is placed in each huse
    
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
    pairs.df[person.index, paste("Group", assigned.group)] = 1
    
    
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
                               progress.bar = TRUE, shuffled.ppl = NA){
    called.fns <<- union(called.fns, "GeneratePairsData")
    if(isTRUE(is.na(shuffled.ppl))){
      stop("No shuffled people")
    }
    col.names = paste("Group", 1:n.groups)
    # pairs.totals = data.frame("Person" = people) %>% 
    #   AddCols(col.names, 0)
    sim.results = vector('list', length = iter.size)
    # loss.results = vector('numeric', length = iter.size)
    # {
    #   nonzero.probs = numeric()
    #   for (i in 1:nrow(prob.df)) {
    #     nonzero.probs[i] = sum(prob.df[i, col.names] != 0)
    #   }
    # } # nonzero probabilities
    # iter.size = round((n.groups + mean(nonzero.probs)) / 2 * sample.size)
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
      # loss.results[i] = loss.fn(sim.pairs, groups, n.groups)
      sim.results[[i]] = AssignPairs(people, groups, n.groups, prob.df, 
                                     nonzero.probs, 1, identical.groups, 
                                     shuffled.ppl = shuffled.ppl)
      if(progress.bar){
        setTxtProgressBar(pb, i)
      }
    }
    if(progress.bar){
      close(pb)
    }
    return(sim.results)
    
    # preserve: nonzero.probs (argument), loss.results (do later)
    # iter.size (replace with sample.size in argumetns)
  }
  OrganizePairs = function(sim.pairs, people, col.names, loss.fn, n.groups = NA){
    called.fns <<- union(called.fns, "OrganizePairs")
    # sim.data = data.frame()
    # for (i in 1:iter.size) {
    #   sim.pairs = sim.results[[i]]
    if(is.na(n.groups)){
      n.groups = length(col.names)
    }
    sim.row = data.frame()
    for (p in 1:length(people)) {
      sim.row[1, paste(people[p], "Group")] = match(1, sim.pairs[p, col.names])
    }
    sim.row$Loss = loss.fn(sim.pairs, groups, n.groups)
    #   sim.data = sim.data %>% rbind(sim.row)
    # }
    # sim.data$Iter = iter.level
    # tryCatch({
    #   all.pairs.data <<- rbind(all.pairs.data, sim.data)
    # }, error = function(e){
    #   all.pairs.data <<- sim.data
    # })
    return(sim.row)
  }
  OrganizeManyPairs = function(sim.results, people, col.names, iter.size, loss.fn, 
                               iter.level = NA, progress.bar = TRUE, n.groups = NA){
    called.fns <<- union(called.fns, "OrganizeManyPairs")
    sim.data = data.frame()
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
    for (i in 1:iter.size) {
      sim.pairs = sim.results[[i]]
      # sim.row = data.frame()
      # for (p in 1:length(people)) {
      #   sim.row[1, paste(people[p], "Group")] = match(1, sim.pairs[p, col.names])
      # }
      # sim.row$Loss = loss.fn(sim.pairs, groups, n.groups)
      sim.row = OrganizePairs(sim.pairs, people, col.names, loss.fn, n.groups)
      sim.data = sim.data %>% rbind(sim.row)
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
    called.fns <<- union(called.fns, "ProbabilityTests")
    # I dont think this depends on the iteration order, 
    # so I am going to leave out shuffled.ppl 
    # print(alpha)
    new.probs = prob.df 
    group.levels = 1:n.groups
    if(identical.groups){
      set.groups = which(prob.df == 1, arr.ind = TRUE)[,2] %>% unique() - 1
      free.groups = setdiff(group.levels, set.groups)
    } else {
      set.groups = group.levels
      free.groups = numeric()
    }
    
    for (p in 1:length(people)) {
      curr.probs = prob.df[p, col.names]
      if(all(curr.probs %in% c(0, 1))){
        new.probs[p, col.names] = curr.probs
      } else {
        group.col.name = paste(people[p], "Group")
        # sim.data[, group.col.name] = sim.data[, group.col.name] %>% 
        #   as.factor()
        # loss.given.group = vector('list', length = length(group.levels))
        # group.levels = levels(sim.data[, group.col.name])
        # unpooled.variance = 0
        # sample.means = numeric()
        
        # mu0 = mean(sim data Loss)
        # xbar = mean(sim data Loss for specific group)
        # h0: mu = mu0 (where mu is derived from xbar)
        # ha: mu < mu0
        
        p.vals = numeric()
        # for (g in group.levels) {
        #   if(g %in% set.groups){
        #     loss.given.group[[g]] = sim.data %>% 
        #       filter_at(group.col.name, function(x) x == g) %>% 
        #       select(Loss) %>% 
        #       unlist()
        #     # unpooled.variance = unpooled.variance + 
        #     #   sd(loss.given.group[[g]]) ^ 2 / length(loss.given.group[[g]])
        #     # sample.means[g] = mean(loss.given.group[[g]])
        #     
        #     # eilminate pairings that associate with worse than average loss values 
        #     if(length(loss.given.group[[g]]) < 1){
        #       p.vals[g] = 0
        #     } else {
        #       tryCatch({
        #         test.result = t.test(loss.given.group[[g]], mu = mean(sim.data$Loss), 
        #                              alternative = 'greater')
        #         p.vals[g] = test.result$p.value
        #       }, error = function(e){
        #         p.vals[g] <<- 0.5
        #       })
        #     }
        #   } else {
        #     p.vals[g] = 0.5
        #   }
        ### delete
        # g6 = numeric()
        # for (p in people) {
        #   g6[p] = mean(sim.data[which(sim.data[, paste(p, "Group")] == 6), "Loss"])
        # }
        # g6 = g6[which(!is.nan(g6))] %>% sort
        ### delete
        for (g in set.groups) {
          p.vals[g] = SingleProbTest(sim.data$Loss, 
                                     which(sim.data[, group.col.name] == g))
        }
        if(length(free.groups) > 0){
          p.vals[(g+1):n.groups] = length(free.groups) * 
            SingleProbTest(sim.data$Loss,
                           which(sim.data[, group.col.name] %in% free.groups))
          # p.vals[(g+1):n.groups] = 0.5
        } # test all free groups as one 
        
        
        # another idea: 
        # hypothesis test for each set group
        # hypothesis test for unset groups as one (all at once)
        # assign to one (or split maybe) set groups IF
        # p.val for unset groups (as a whole) is significant
        # p.val for (at least one) set group is significant (for less than)
        
        # if(p.vals[[g]] > pValueThreshold){
        #   p.vals[[g]] = 1
        # }
        # }
        # new.probs[p, col.names] = (1 - p.vals) / sum(1 - p.vals)
        new.probs[p, col.names] = ifelse(p.vals < alpha, 0, 1) %>% 
          ProbVector()
      }
    }
    return(new.probs)
  }
  SingleProbTest = function(loss.data, test.index){
    called.fns <<- union(called.fns, "SingleProbTest")
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
    called.fns <<- union(called.fns, "CalcNonzeroProbs")
    if(isTRUE(is.na(col.names))){
      col.names = names(prob.df[-1])
    }
    nonzero.probs = numeric()
    for (i in 1:nrow(prob.df)) {
      nonzero.probs[i] = sum(prob.df[i, col.names] != 0)
    }
    return(nonzero.probs)
  }
  FindOverfullGroups = function(prob.df, groups, n.groups){
    called.fns <<- union(called.fns, "FindOverfullGroups")
    overfull.groups = numeric()
    for (i in 1:n.groups) {
      if(sum(prob.df[, paste("Group", i)] == 1) > groups[i]){
        overfull.groups = c(overfull.groups, i)
      }
    }
    return(overfull.groups)
  }
} # functions
folder.path = "/Users/malexk999/Desktop/Cloud desktop/University of Maryland/Frat/Big little 2021/Algorithm explanation"
file.name = "Example data.xlsx"
twins.elig = c()
# twins.elig = "all"
# twins.elig = c("Howard", "Samuel")

# make it so that littles MUST be paired with a big they requested
use.default.groups.shortcut = FALSE # TRUE or FALSE
# 1 minute for TRUE, 2 minutes for FALSE

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
  
  loss.fn = function(pairs.df, groups, n.groups){
    total.loss = 0
    for (i in 1:n.groups) {
      group.members = pairs.df %>% 
        filter_at(paste("Group", i), function(x) x > 0) %>% 
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
    answer = MakePairs(people, groups, loss.fn, default.groups = default.groups)
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
    # loss.fn(example.pairing, groups, ncol(example.pairing) - 1) # 66
    
    
  } # example 
} # where the magic happens 
