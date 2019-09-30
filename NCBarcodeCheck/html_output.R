## Html output of experiment_design



library('htmlTable')
library("dplyr")

# The styles and format to apply to the result table ----------------------

# nucleotide
style = "margin-bottom : 50px;
margin-top: 1000px;color : black;
border: solid 2px black;
border-collapse: collapse; 
heigt : 30px; width : 50px;
background-color:"

green = paste (style,"#99e7a9")
yellow = paste (style,"#ffcf40")
red = paste(style, "#f16a6a")
blue = paste(style,"#6BB8E7")
blackish = paste(style, "#999999")

# sample column or Id column
general_style = "text-align:center;
background-color: white;
color : black; 
font-family : verdana; font-size : 90%;
padding: 5px 20px 5px 20px;"

cgroup_style = "text-align:center;
background-color: white;
color : black;
border-top : 2px solid grey;
border-bottom : 2px solid white;
padding : 5px 5px 20px 5px;
font-family : verdana;
border-collapse : collapse;"

css.cgroup = c(cgroup_style, cgroup_style, cgroup_style)
css.rgroup.sep = "border-bottom : 2px solid black;"
css.rgroup = "border-bottom: 2px solid black; 
font-family : verdana;
font-size ; 60%;
font-style : oblique;
padding-top : 20px;"


# The functions -----------------------------------------------------------



# styles according to platform
table_style = function(sequence_m, platform){
  if (platform == 0){
    style_sequence = sequence_m %>% 
      gsub("A",green, .)%>%
      gsub("C", blue, .)%>%
      gsub("T", red, .)%>%
      gsub("G", blackish, .) 
  }else if (platform == 1){
    style_sequence = sequence_m %>% 
      gsub("A",green, .)%>%
      gsub("C", red, .)%>%
      gsub("T", yellow, .)%>%
      gsub("G", blackish, .) 
  }else if (platform == 2){
    style_sequence = sequence_m %>% 
      gsub("A",yellow, .)%>%
      gsub("C", red, .)%>%
      gsub("T", green, .)%>%
      gsub("G", blackish, .) 
  }else if (platform == 4){
    style_sequence = sequence_m %>% 
      gsub("A|C", red, .)%>%
      gsub("G|T", green, .)
  }
  return(style_sequence)
}





# convert a df column contening sequence into nucleotide matrix
split_sequence = function(df_column){
  col = nchar(df_column[1])
  row =length(df_column)
  column = as.matrix( unlist(strsplit(df_column,"")))
  column = matrix(column, nrow = row, ncol = col, byrow = TRUE)
  return (column)
}


# The final output according to the type of platform
build_table_style = function(result_df, platform){
  if(ncol(result_df) != 4 || 6){
    if(ncol(result_df) == 4){
      return(build_table_style_single(result_df, platform))
    }else if (ncol(result_df == 6)){
      return(build_table_style_dual(result_df, platform))
    } 
  }else{
    return (NULL)
  }
}


# html_output for dual indexing
build_table_style_dual = function(result_df, platform){
  nb_lane = as.numeric(result_df$Lane[nrow(result_df)])
  multiplexing_level = nrow(result_df)/nb_lane
  rgroup = paste("Lane", as.character(c(1: nb_lane)))
  n.rgroup = c(rep(multiplexing_level,nb_lane))
  cgroup = c("Id1", "Sequence 1", "Sample", "Sequence 2", "Id2")
  css.cgroup = rep(cgroup_style,5)
  result_m = as.matrix(result_df)
  nucleotide_m1 = split_sequence(result_df$sequence1)
  nucleotide_m2 = split_sequence(result_df$sequence2)
  general_style_m = matrix(general_style,
                           nrow= nrow(result_m), ncol = 1) 
  css.cell = cbind(general_style_m,
                   table_style(nucleotide_m1,platform),
                   general_style_m,
                   table_style(nucleotide_m2,platform),
                   general_style_m)
  result_m = cbind(result_m[,"Id1"],
                   nucleotide_m1, result_m[,"sample"],
                   nucleotide_m2, result_m[,"Id2"])
  n.cgroup = c(1, ncol(nucleotide_m1), 1,
               ncol(nucleotide_m2), 1)
  result_m = htmlTable(result_m, 
                       cgroup = cgroup,
                       n.cgroup = n.cgroup,
                       css.cgroup =css.cgroup,
                       rgroup = rgroup,
                       n.rgroup = n.rgroup,
                       css.rgroup.sep = css.rgroup.sep,
                       css.rgroup = css.rgroup,
                       css.cell = css.cell
  )
  return (result_m)
}


# html_output for single indexing
build_table_style_single = function(result_df, platform){
  nb_lane = as.numeric(result_df$Lane[nrow(result_df)])
  multiplexing_level = nrow(result_df)/nb_lane
  rgroup = paste("Lane", as.character(c(1: nb_lane)))
  n.rgroup = c(rep(multiplexing_level,nb_lane))
  cgroup = c("Sample", "Sequence", "Id")
  result_m = as.matrix(result_df)
  nucleotide_m = split_sequence(result_df$sequence)
  general_style_m = matrix(general_style,
                           nrow= nrow(result_m),
                           ncol = 1) 
  css.cell = cbind(general_style_m,
                   table_style(nucleotide_m,platform),
                   general_style_m)
  result_m = cbind(result_m[,"sample"], nucleotide_m, result_m[,"Id"])
  n.cgroup = c(1, ncol(nucleotide_m), 1)
  result_m = htmlTable(result_m, 
                       cgroup = cgroup,
                       n.cgroup = n.cgroup,
                       css.cgroup =css.cgroup,
                       rgroup = rgroup,
                       n.rgroup = n.rgroup,
                       css.rgroup.sep = css.rgroup.sep,
                       css.rgroup = css.rgroup,
                       css.cell = css.cell
  )
  return (result_m)
}

