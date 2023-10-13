linear_predictor_formula <- function(formula = 'base',            #Which formula
                                     temporal_model = 'bym',      #What model does temporal effecst follow
                                     temporal_rank_def = 1,       #How rank deficient is the temporal effects
                                     temporal_structure_matrix,   #The temporal effects structure matrix
                                     temporal_hyper,              #Hyper parameters and priors for temporal effects
                                     spatial_model = 'bym',       #What model does spatial effects follow
                                     spatial_rank_def = 1,        #How rank deficient is spatial effects
                                     spatial_structure_matrix,    #The spatial effects structure matrix
                                     spatial_hyper){              #Spatial effects hyper parameters and hyper priors
  
  
  #Add different models based on structure, besag vs bym2, etc...
  allowed_formulas = c('pure_iid', 'base', 'typeI', 'typeII', 'typeIII', 'typeIV')
  
  #Check that specified formula is within allowed_formulas, then continue
  if(formula %in% allowed_formulas){
    
    #Make the base formula (is the overall mean specified correctly here???)
    base_formula <- deaths ~ 1 + f(year, 
                                   model = temporal_model,
                                   scale.model = T, 
                                   constr = T, 
                                   rankdef = temporal_rank_def,
                                   graph = temporal_structure_matrix,
                                   hyper = temporal_hyper
                                   ) + 
                                  f(county, 
                                    model = spatial_model,
                                    scale.model = T,
                                    constr = T,
                                    rankdef = spatial_rank_def,
                                    graph = spatial_structure_matrix,
                                    hyper = spatial_hyper
                                  )
      
    
    
    if(formula == 'base'){
      #Return the base formula
      return(base_formula)
      
    } else if(formula == 'typeI'){
      #Add typeI interaction to base fomrula
      
      #Return formula w. typeI interaction
      
    } else if(formula == 'typeII'){
      #Add typeII interaction to base fomrula
      
      #Return formula w. typeII interaction
    } else if(formula == 'typeIII'){
      #Add typeIII interaction to base fomrula
      
      #Return formula w. typeIII interaction
    } else if(formula == 'typeIV'){
      #Add typeIV interaction to base fomrula
      
      #Return formula w. typeIV interaction
    }
    
    
  } else{ #Otherwise print that model not allowed
    print("Not a permissible formula")
  }
}
