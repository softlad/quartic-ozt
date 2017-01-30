package com.matpollard.speech;

public Class Utils {
	
	public static int whatSection(Section[] sections, int time)	{
	  
	  if (time<0)
	  {
	    System.out.println("error in What Section");
	    exit(0);
	  }
	  
	  int section=-1;
	  
	  for (int i=0; i<nsections.length; i++) {
		 if ((section[i].beg<=time) && (section[i].end>=time)) {
		 	section=i;
		 }
	  }
	    
	  return (section);
	}

	public static int getBlockLength(int time, Sections[] sections) { 
  
  		short avperiod;
  
  		int section = whatSection(sections, time);
  
  		if (section != -1) {
  			avperiod=section[section].avperiod;
  		}
  		else {
  			int sum = 0;
  			for (int s=0; s < sections.length; s++) {
  				sum += section[s].avperiod;
  			}
  			avperiod=round((double)sum/sections.length);
  		}
  
  		int length = (2.5 * avperiod);

  		if (length%2==0) {
  			length++;
  		} 
	}


}