				 
var g_numCheckedCheckBoxes=0;

function SetFormSubmitEvent()
{
	jQuery("#alignments").find("input[name='RESUB_SEQLOC']").each(function (index) {        		
		jQuery(this).bind("click",initSubmButton);        
	});	
	jQuery("#submReanl2").bind("click", function (e) {jQuery("#submReanl1").click();});
	jQuery("#reanalysisform").bind("submit",reAnalyze);
}

function reAnalyze(e) {
	var intendTargets = "";
	jQuery("#alignments").find("input[name='RESUB_SEQLOC']").each(function (index) {        		
		if(this.checked) {
			if(intendTargets == "") {
				intendTargets += this.value;
			}
			else if (intendTargets.indexOf(this.value) == -1) {
				intendTargets += "," + this.value;
			}
		}		
	});	
	var submitTargets = intendTargets.split(",");
	for(i = 0;i < submitTargets.length;i++) {
        elem = document.createElement("input");
        jQuery(elem).attr("type", "hidden");
        jQuery(elem).attr("name", "RESUB_SEQLOC");
        jQuery(elem).val(submitTargets[i]);                    
        jQuery("#reanalysisform").append(elem);   
        //console.log("APPEND: " + elem.name + " value=" + elem.value);                         
    }    	      
}

function CheckTheSameTarget(originalTarget)
{
	jQuery("#alignments").find("input[name='RESUB_SEQLOC']").each(function (index) {        		
		if(this.value == originalTarget.value) {			
			this.checked = originalTarget.checked;
		}
	});	
	return false;

}
function initSubmButton(e)
{
	if(this.checked) {		
		g_numCheckedCheckBoxes++;		
		jQuery("#submReanl1").removeAttr("disabled");
		jQuery("#submReanl2").removeAttr("disabled");		
	}
	else {
		g_numCheckedCheckBoxes--;
		if(g_numCheckedCheckBoxes == 0) {
			jQuery("#submReanl1").attr("disabled","disabled");
			jQuery("#submReanl2").attr("disabled","disabled");			
		}		
	}
	CheckTheSameTarget(this);	
}

utils.addEvent(window,
				"load", SetFormSubmitEvent,				
				false);


