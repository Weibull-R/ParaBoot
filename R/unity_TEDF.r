# this function is required to fill in censored event times in the unity_TEDF					
unity_TEDF<- function(pivotal_pts, event_vec) {					
	time_vec<-NULL				
	piv_num<-1				
					
	for(ev in 1:length(event_vec)) {				
		if(event_vec[ev] == 1) {			
			time_vec <- c(time_vec, pivotal_pts[piv_num])		
			piv_num <- piv_num + 1		
		}else{			
			if(piv_num > 1 && piv_num<length(pivotal_pts)+1) {		
				midtime<-sum(pivotal_pts[(piv_num-1):piv_num])/2	
				time_vec <- c(time_vec, midtime)	
			}		
			if(piv_num == 1 ) {		
				midtime<- pivotal_pts[1]/2	
				time_vec <- c(time_vec, midtime)	
			}		
			if(piv_num >length(pivotal_pts) )  {		
				endtime<-max(pivotal_pts) + 10	
				time_vec <- c(time_vec, endtime)	
			}		
		}			
	}				
	return(data.frame(time=time_vec, event=event_vec))				
}