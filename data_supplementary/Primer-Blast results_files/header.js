

        // add a back_url to a#myncbi's href
        /*
        var acctLinks = document.querySelectorAll('#account_login, #account_logout');
        for (var i = 0; i < acctLinks.length; i++) {
            acctLinks[i].href = acctLinks[i].href + '?back_url=' + document.location.href;
        }
        */
        var cubby = getCookie('WebCubbyUser');
        cubby = decodeURIComponent(decodeURIComponent(cubby));

        var username = getUser(cubby);

        if (username) {
            jQuery('#uname_short').text(username.trunc(20));
            jQuery('#uname_long').text(username.trunc(40));
            jQuery('#account_login').hide();
            jQuery('#account_info').show();
        }
        else {
            jQuery('#account_login').show();
            jQuery('#account_info').hide();
        }


        function getUser(c) {

            var re_logd = /.*logged-in\=(\w*);.*/;
            var re_user = /.*my-name\=([\w|\-|\.|\ |\@|\+]*);.*/;
            if (c) {
                var l = re_logd.exec(c);
                if (l && l[1] && l[1] === 'true') {
                    var u = re_user.exec(c);
                    if (u && u[1]) {
                        return u[1];
                    }
                }
            }
            return '';
        }


        function getCookie(f) {
            var e;
            if (window.sessionStorage) {
                try {
                    e = sessionStorage.getItem(f) || '';
                } catch (g) {
                    e = '';
                }
                if (e.length > 0) {
                    return e;
                }
            }
            if (document.cookie.length > 0) {
                e = document.cookie.indexOf(f + '=');
                if (e !== -1) {
                    e = e + f.length + 1;
                    f = document.cookie.indexOf(';', e);
                    if (f === -1) {
                        f = document.cookie.length;
                    }
                    return unescape(document.cookie.substring(e, f))
                }
            }
            return '';
        }


        function updateUserSurvey(userSurveyName) {

            userSurveyUrl = location.protocol + "//" + location.host + "/blast/fmtsave.cgi";
            var rp = new RemoteDataProvider(userSurveyUrl);

            rp.onSuccess = function (obj) {
                console.log("updateUserSurvey: fmtsave.cgi success, params: " + params);
            };
            rp.onError = function(obj) {
                console.log("updateUserSurvey: fmtsave.cgi error, params: " + params);
            }

            var params = "CMD=Update&CUBBY_OBJECT=UserSurvey&USER_SURVEY_NAME=" + userSurveyName + "&BANNER_AD_WAS_SHOWN=true"
            rp.Request(params);
        }

        function RetrieveUserSurveyInfo(userSurveyName) {
            userSurveyUrl = location.protocol + "//" + location.host + "/blast/fmtsave.cgi";

            var rp = new RemoteDataProvider(userSurveyUrl);

            rp.onSuccess = function(obj) {
                if(obj.responseText == "0") {
                    jQuery("#banner_ad").show();
                }
            };
            rp.onError = function(obj) {
                console.log("fmtsave.cgi: Error saving " + this.name + " configuration");
            }

            var params = "CMD=Get&CUBBY_OBJECT=UserSurvey&USER_SURVEY_NAME=" + userSurveyName;
            rp.Request(params);
        }

function bannerAdInit()
{
    var userSurveyName = "";
    //For now we only specify "blastspec", pageType/program may be used in future
    var pagetype = jQuery("#closeAd").attr("pagetype");
    var program = jQuery("#closeAd").attr("program");
    var blastSpec = jQuery("#closeAd").attr("blastspec");
    //Construct survey name <page_type>_<program>_<blast_spec>
    //For example: BlastSearch_blastp_blast2seq or BlastResults_blastn
    if(pagetype) {
        userSurveyName = pagetype;
    }
    if(program) {
        userSurveyName += (userSurveyName == "") ? program : "_" + program;
    }
    if(blastSpec) {
        userSurveyName += (userSurveyName == "") ? blastSpec : "_" + blastSpec;
    }
    jQuery("#banner_ad #closeAd").bind("click", function (e) {
        var ariaContrl = jQuery(this).attr("aria-controls");
        jQuery("#" + ariaContrl).remove();
        updateUserSurvey(userSurveyName);
    });
    jQuery("#banner_ad #jt").bind("click", function (e) {
        jQuery("#btnJt").click();
    });
    jQuery("#btnJt").bind("click", function(e) {
        winRef = window.open(jQuery("#jt")[0].href, "hlp" )
        jQuery("#closeAd").click();
    });
    RetrieveUserSurveyInfo(userSurveyName);
}
        var $popupMenu = jQuery('#account_info').ncbipopup();
        if(jQuery("#banner_ad")[0]) bannerAdInit();


