function(in_str){
  if(is.na(in_str)){return(3)} else if(strcmp(in_str,"Positive")){
    return(1) } else if(strcmp(in_str,"Negative")){
      return(2) }else if(strcmp(in_str,"Missing")){
        return(3) }else {
          return(3) }
  
}
