
/*  $Id: loader.js 47441 2023-04-05 17:00:08Z rudnev $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * File Description: Dynamic scripts loader and synchronizer
 */

timeStamp = new Date().getTime();
(function(){
    var isIE = (window.navigator.userAgent.indexOf('Trident/') >= 0);
    var isGooglebot = (window.navigator.userAgent.indexOf('Googlebot/') >= 0);
    var refNode = document.currentScript;
/*    var refNode = document.scripts[document.scripts.length - 1];
    if (!refNode || !refNode.src) { // SV-4715 "Standalone SV broken (Brave)"
        window.alert('Not compatible browser!');
        return;
    }*/
    var srcjs = refNode.src.split('/');
    var viewer = ({sviewer: 'SeqView', treeviewer: 'TreeView', multialign: 'MultiAlignView'})[srcjs.pop().replace('.js', '')];
    if (window[viewer]) return;
    var dl = document.location;
    var externJS = (srcjs[2] != dl.hostname ? '-extern.js' : '.js');
    var dotJS = (srcjs[2].search(/www|qa/) !== 0 || dl.href.indexOf('debug') > 0) ? '-debug.js' : '.js';
    
    var domNCBI = 'ncbi.nlm.nih.gov';
    var re_https_only_site = /(www|dev|test)\.ncbi\.nlm\.nih\.gov/;
    var webNCBI = (dl.hostname.substr(-domNCBI.length) == domNCBI &&
                   !re_https_only_site.test(srcjs[2]) ? dl.protocol : 'https:') + '//'
        + (srcjs[2].indexOf(domNCBI) == -1 || srcjs[2].indexOf('blast.' + domNCBI) >= 0
        ? ('www.' + domNCBI) : srcjs[2]) + '/';

    var ExtJSver = dl.href.substr(dl.href.search(/extjs/), 10).slice(5) || '7.1.0';
    if (!document.querySelector('meta[name=viewport]')) {
        var meta = document.createElement('meta');
        meta.setAttribute('name', 'viewport');
        meta.setAttribute('content', 'width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no');
        document.getElementsByTagName('head')[0].appendChild(meta);
    }

    var appPath = ''
    for (var i = 3; i < srcjs.length - 1; i++)
        appPath += srcjs[i] + (srcjs[i] ? '/' : '');
    var extPath = webNCBI + 'core/extjs/ext-' + ExtJSver + '/build/';
//    if (dl.href.search(/extjs7./) > 0) extPath = 'https://dev.ncbi.nlm.nih.gov/staff/borodine/ext-' + ExtJSver + '/build/';
    var gbPath = appPath;
    var thirdPartyScripts = {
        Ext: {tag: 'script', attr: {type: 'text/javascript', src: extPath + 'ext-all' + dotJS}},
        jQuery: {tag: 'script', attr: {type: 'text/javascript', src: webNCBI + 'core/jquery/jquery-3.6.0.js'}},
        UUD: {tag: 'script',  attr: {type: 'text/javascript', src: webNCBI +'projects/genome/uud/js/uud' + externJS}},
        TMS: {tag: 'script', attr: {type: 'text/javascript', src: webNCBI + 'projects/genome/trackmgr/0.7/js/tms' + externJS}},
        GB: {tag: 'script', attr: {type: 'text/javascript', src: webNCBI + gbPath + 'js/gbx-hub' + (isIE ? '.js' : dotJS)}},
        ncbi: {tag: 'script', attr: {type: 'text/javascript', src: webNCBI + 'portal/portal3rc.fcgi/rlib/js/InstrumentOmnitureBaseJS/InstrumentNCBIBaseJS/InstrumentPageStarterJS.js'}},
        MSAwasm: {tag: 'script', attr: {type: 'text/javascript', src: webNCBI + appPath + 'js/align_init.js', async: true}},
        ExtCSS: {tag: 'link', attr: {rel: 'stylesheet', type: 'text/css', href: extPath + 'classic/theme-gray/resources/theme-gray-all.css'}},
        GBCSS: {tag: 'link', attr: {rel: 'stylesheet', type: 'text/css', href: webNCBI + gbPath + 'css/gbc/gb-hub.css'}}, 
        CSS: {tag: 'link', attr: {rel: 'stylesheet', type: 'text/css', href: webNCBI + appPath + 'css/style.css'}}};

    window[viewer + 'OnReady'] = function(callback, scope) {
        if (typeof window[viewer] === 'undefined')
            setTimeout(function() { window[viewer + 'OnReady'](callback, scope); }, 100);
        else {
            Ext.ariaWarn = Ext.emptyFn;
            Ext.onReady(callback, scope);
        }
    }

    var insertElement = function(resName, callback) {
        var res = thirdPartyScripts[resName];
        var deferCallback = !callback ? function(){} : function() {
            if (dotJS.length > 3) console.log(resName);
            if (!window[resName]) setTimeout(deferCallback, 30);
            else callback();
        }

        if (res.tag == 'link') {
            for (var p in document.styleSheets)
                if (res.attr.href == document.styleSheets[p].href) { res = false; break; }
        } else {
            if (typeof window[resName] !== 'undefined') res = false;
            else {
                for (var p in document.scripts)
                    if (res.attr.src == document.scripts[p].src) { res = false; break;}
            }
        }
        if (!res) {
           deferCallback();
           return;
        }

        var el = document.createElement(res.tag);
        for (var p in res.attr) el[p] = res.attr[p];
        if (typeof callback === 'function') {
            el.onload = el.onreadystatechange = function() {
                if (!this.readyState || this.readyState === "loaded" || this.readyState === "complete") {
                    if (isIE)
                        try {
                            var tmp = document.body.appendChild(document.createElement('div'));
                            tmp.innerText = ' ';
                            document.body.removeChild(tmp);
                        } catch(e){};
                    this.onreadystatechange = this.onload = null;
                    deferCallback();
                }
            }
        }
        refNode.parentNode.insertBefore(el, refNode);
    }
    insertElement('ExtCSS');
    insertElement('CSS');
    if (refNode.getAttribute('extjs') == 'skip') window.Ext = window.Ext || {};

    var counter = 100;

    insertElement('Ext', function finalLoad() {
        if ((typeof Ext === 'undefined' || !Ext.onReady) && counter--) {
            setTimeout(finalLoad, 100);
            return;
        }
        if (Ext.getVersion && Ext.getVersion().major >= 5) {
            Ext.enableAriaButtons = false;
            Ext.scopeCss = true;
            Ext.namespace(viewer);
            Ext.apply(window[viewer], {webNCBI: webNCBI, base_url: webNCBI + appPath, standalone: dl.pathname.indexOf(appPath) == 1});
            window['init' + viewer]();
        } else {
            alert('Current version of the Application works only in ExtJS ver. 6.+ environment!');
            return;
        }
        if (Ext.isIE9) thirdPartyScripts['UUD'].attr.src = thirdPartyScripts['UUD'].attr.src.replace(externJS, '-extern.js');
        if (document.location.hostname.indexOf(domNCBI) == -1) {
            window[viewer].origHostname = document.location.hostname;
            __ncbi_stat_url = 'https://www.ncbi.nlm.nih.gov/stat';
            ncbi_signinUrl = 'https://www.ncbi.nlm.nih.gov/portal/signin.fcgi?cmd=Nop';
            ncbi_pingWithImage = true;    
        } else {
            if (webNCBI.indexOf('https://') != 0) {
               __ncbi_stat_url = 'https://dev.ncbi.nlm.nih.gov/stat';
               ncbi_pingWithImage = true;
            }
        }
        if (viewer == 'MultiAlignView' && !Ext.isIE && !isGooglebot) {
//            DataReady = function() { console.log(arguments, this);}
            insertElement('MSAwasm');
        }
        insertElement('jQuery', function(){
            insertElement('UUD', function() { if (viewer == 'SeqView') {
                insertElement('TMS');
                insertElement('GB');
                insertElement('GBCSS');
            }});
        });
        if (!Ext.removeNodeOrig) {
            Ext.removeNodeOrig = Ext.removeNode;
            Ext.removeNode = function(n) { if (n) Ext.removeNodeOrig(n); }
        }
        if (!(Ext.isEdge|Ext.isSafari|Ext.isFirefox|Ext.isChrome|Ext.isChromeiOS|Ext.isChromeMobile|isGooglebot)) {
            alert('Warning!\nUnsupported internet browser. Some functions may not work.');    
        }

        window[viewer].extPath = extPath;
        window[viewer].extImagesPath = extPath + 'classic/theme-gray/resources/images/';
        window[viewer + 'OnReady'](function(){
            insertElement('ncbi');
            window[viewer].loadApp(refNode);
        });
    });
    
    if (window.NCBIGBUtils && NCBIGBUtils.makeTinyURL) return;

    NCBIGBUtils = window.NCBIGBUtils || {};
    NCBIGBUtils.makeTinyURL = function(url, callback) {
        callback('');
/*        if (url.indexOf('http') != 0) {
            callback(url);
            return;
        }
        var key = '9add369550da0ca4a4e899358f70be53';
        var url_to_pass = 'https://go.usa.gov/api/shorten.json?login=NLM%20NCBI&apiKey=' + key + '&longUrl=' +
            encodeURIComponent(url);

        jQuery.get(url_to_pass, function(data, status, xhr) {
            var res = 'go.usa.gov url shortener ';
            try {
                res = (data.response.data.entry[0].short_url || '').replace(/http:/, 'https:');
            } catch (e) {
                if (data && data.response && data.response.errorMessage)
                    res += ': ' + data.response.errorMessage;
                else
                    res += 'is not avallable at this moment';
                if (ncbi && ncbi.sg) ncbi.sg.ping({ event: 'url_shortener_failure', err: res });
            }
            callback(res);
        });*/
    }

})();




/*  $Id: globals.js 46457 2021-05-24 13:25:18Z borodine $
 * Authors:  Maxim Didenko
 * File Description:
 *
 */

function initSeqView()
{

SeqView.doInitPing = true;
SeqView.loadApp = function (refNode) {
    var items = Ext.query('div[class=SeqViewerApp]');
    for (var i = 0; i < items.length; i++) {
        var item = items[i];
        var id = item.id || Ext.id();
        if (!item.hasAttribute('data-autoload') && refNode.id != 'autoload') continue; 
        if (SeqView.App.findAppByDivId(id)) return; 
        item.id = id;
        var app = new SeqView.App(id);
        app.load();
    }
}

if (Ext.isIE) {
    var db = document.querySelector('meta[name=ncbi_db]');
    SeqView.geneIE = (db && db.content && (db.content == 'clone' || db.content == 'gene'));
}

Ext.define('SeqView.SearchField', {
    extend: 'Ext.form.field.Text',
    width:180,
    triggers: { 
        clear: {
            cls: 'x-form-clear-trigger',
            hidden: true,
            handler: function(tb, trg) { this.clickTrigger(trg.id); }
        },
        search: {
            cls: 'x-form-search-trigger',
            handler: function(tb, trg) { this.clickTrigger(trg.id); }
        }
    },
    initComponent: function(){
        this.callParent(arguments);
        this.on('specialkey', function(f, e){
            if(e.getKey() == e.ENTER){
                this.clickTrigger(null, {id: 'search'});
            }
        }, this);
    },

    clickTrigger: function(trgid) {
        if (this.callback) this.callback();
        if (trgid == 'clear') {
            this.triggers.clear.hide();
            this.setValue('');
        } else this.triggers.clear.show();
        if (this.userClickTrigger) {
            this.userClickTrigger(this.getRawValue());
            return;
        }
        var bp = this.store.baseParams;
        bp = bp || {};
        bp.query = this.getRawValue();
        bp.regexp = bp.query ? globStringToRegex(bp.query) : '';
        if (bp.filterBy) {
            this.store.filterBy(bp.filterBy, this.store);
        } else this.store.reload({params:{start: 0}});
        
    }
});

if (navigator.maxTouchPoints > 1 || 'ontouchstart' in document) {
    SeqView.useTouch = true;
}

SeqView.Cookies = function(config) {
    this.path = "/";
    this.expires = new Date(new Date().getTime()+(1000*60*60*24*7)); //7 days
    this.domain = ".nih.gov"; // this is essential for proper cookie functioning
    this.secure = false;
    Ext.apply(this, config);

    var readCookies = function() {
        var cookies = {};
        var c = document.cookie + ";";
        var re = /\s?(.*?)=(.*?);/g;
        var matches;
        while((matches = re.exec(c)) != null){
            var name = matches[1];
            // fix incorect cookie name for IE8 (like 'WebCubbyUser; sv-user-data')
            var a = name.split(';');
            name = a[a.length-1].trim();
            var value = matches[2];
            cookies[name] = value;
        }
        return cookies;
    }
    this.state = readCookies();
};

SeqView.Cookies.MaxDate = new Date(new Date('10000/01/01 GMT').getTime()-100);
SeqView.Cookies.UserTracksCookieNameBase = 'sv-usertracks-key';
SeqView.Cookies.UserDataCookieNameBase = 'sv-userdata-key';
SeqView.Cookies.AppDataCookieNameBase = 'sv-data-key';
SeqView.Cookies.UserTracksCookieName = SeqView.Cookies.UserTracksCookieNameBase;
SeqView.Cookies.UserDataCookieName = SeqView.Cookies.UserDataCookieNameBase;
SeqView.Cookies.AppDataCookieName = SeqView.Cookies.AppDataCookieNameBase;

SeqView.Cookies.prototype = {
    get : function(name, defaultValue){
        return typeof this.state[name] == "undefined" ?
            defaultValue : this.state[name];
    },

    set : function(name, value){
        if(typeof value == "undefined" || value === null){
            this.clear(name);
            return;
        }
        this.setCookie(name, value);
        this.state[name] = value;
    },

    clear : function(name){
        delete this.state[name];
        this.clearCookie(name);
    },

    // private
    setCookie : function(name, value){
        document.cookie = name + "=" + value +
           ((this.expires == null) ? "" : ("; expires=" + this.expires.toGMTString())) +
           ((this.path == null) ? "" : ("; path=" + this.path)) +
           ((this.domain == null) ? "" : ("; domain=" + this.domain)) +
           ((this.secure == true) ? "; secure" : "");
    },

    // private
    clearCookie : function(name){
        document.cookie = name + "=del; expires=Thu, 01-Jan-1970 00:00:01 GMT" +
           ((this.path == null) ? "" : ("; path=" + this.path)) +
           ((this.domain == null) ? "" : ("; domain=" + this.domain)) +
           ((this.secure == true) ? "; secure" : "");
    }
};
SeqView.SessionData = new SeqView.Cookies({expires:null, secure:Ext.isSecure});
SeqView.UserTracks = new SeqView.Cookies({expires:SeqView.Cookies.MaxDate, secure:Ext.isSecure});
SeqView.trackSets = {tmsSets: [], usrSets: [], trackLists: []};

SeqView.configureTrackSets = function(app) {
    TMS.TrackSets.TrackSetService.SetAssembly(app.m_AssmContext);
    TMS.TrackSets.GetDefaultTrackSet('SViewer', app.m_AppContext, app.m_AssmContext)
        .done(function() {
            app.m_defaultTrackSet = this.GetTracks();
        })
        .fail(function() {
            console.log('Failed to get default track set');
        });
    SeqView.requestTrackSets = function(callback) {
        TMS.TrackSets.TrackSetService.GetTracksets(true)
           .done(callback)
           .fail(function() { callback([]);}
        );
    }
};

SeqView.AreaFlags = {
    Link:           (1 << 0),  // a direct link stored in m_Action
    CheckBox:       (1 << 1),  // a toggle button/icon
    NoSelection:    (1 << 2),  // the object can't be selected
    ClientSelection:(1 << 3),  // the selection can be done on client
    NoHighlight:    (1 << 4),  // on highlighting on mouse over
    NoTooltip:      (1 << 5),  // do not request and show tooltip
    TooltipEmbedded:(1 << 6),  // tooltip embedded
    Track:          (1 << 7),  // track title bar
    Ruler:          (1 << 8),  // ruler bar
    Editable:       (1 << 9),  // editable area
    NoPin:          (1 << 10), // not pinnable
    IgnoreConflict: (1 << 11), // feature can be ignored (isca browser feature editing only)
    Sequence:       (1 << 12), // sequence track flag
    Comment:        (1 << 13), // render a label/comment on client side
    DrawBackground: (1 << 14), // highlight background for this area
    Dirty:          (1 << 16), // dirty flag
    NoNavigation:   (1 << 17), // no havigation buttons on title bar
    Legend:         (1 << 18), // legends for graph_overlay tracks
    NoCaching:      (1 << 19), // The tooltip for this feature should not be cached
    ZoomInOnDblClck:(1 << 20)  // Zoom in tracks on double click event
};


SeqView.getHelpURL = function(add) {
    return 'https://www.ncbi.nlm.nih.gov/tools/sviewer/' + (add || '');  
};

SeqView.showHelpDlg = function(extra_params) {
    window.open(SeqView.getHelpURL() + (extra_params ? extra_params : ''));
};

SeqView.showAboutMessage = function() {
    var msg = '<p>NCBI Graphical Sequence Viewer - graphical display for the Nucleotide and Protein sequences.</p>Version '
    + SeqView.getVersionString() + '<br><br><span style="font-size:10px;">' + SeqView.version.rev + '</span><br>';
    msg += '<br>ExtJS version: ' + Ext.getVersion('core').version + '<br>';
    msg += '<br><a href=\"' + SeqView.base_url + 'info.html' + '\" target=\"_blank\" style=\"color:blue\">CGI binaries Info</a>';
    msg += '<br><br><a href=\"'+' https://www.ncbi.nlm.nih.gov/tools/sviewer/' + '\" target=\"_blank\" style=\"color:blue\">Documentation</a>';
    msg += '<br><br><a href=\"'+' https://www.ncbi.nlm.nih.gov/tools/gbench/' + '\" target=\"_blank\" style=\"color:blue\">Graphical Sequence Viewer is also available as part of NCBI Genome Workbench</a>';
    Ext.Msg.show({
        title: 'Sequence Viewer',
        msg: msg,
        maxWidth: '350',
        minWidth: '300',
        buttons: Ext.Msg.OK
    });
}

SeqView.getVersion = function() { return parseFloat(SeqView.version.num); }
SeqView.getVersionString = function() { return SeqView.version.str }

SeqView.decode = function(data) { return (typeof data === 'object') ? data : Ext.decode(data); }

SeqView.m_Apps = [];

SeqView.tmsSearch = async function(params, callback) {
    const data = await TMS.Direct.search(params);
    callback(data);
} 

SeqView.version = {rev: 'Revision: 47578, install date: 2023-10-05 14:38', str: '3.49.0', num:'3490'};
/* $Id: utils.js 45694 2020-10-26 15:17:12Z borodine $
 * File Description: Univirsal tools and objects
 */
if (!window.NCBIGBUtils || !window.NCBIGBUtils.ClearBrowserSelection) {
var utils = window.NCBIGBUtils = window.NCBIGBUtils || {};
// Decorate long numbers with commas
Number.prototype.commify = function() {
    nStr = this + '';
    x = nStr.split('.');
    x1 = x[0];
    x2 = x.length > 1 ? '.' + x[1] : '';
    var rgx = /(\d+)(\d{3})/;
    while (rgx.test(x1)) {
        x1 = x1.replace(rgx, '$1' + ',' + '$2');
    }
    return x1 + x2;
};

Number.prototype.shorten = function() {
    var value = this.valueOf();
    var negv = ( value < 0 );
    if( negv ) value = -value;
    var Suffixes = [ '', 'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y' ];
    var sfx = '';
    for( var i = 0; i < Suffixes.length; i++ ){
        if( value < 1000 ){
            sfx = Suffixes[i];
            break;
        }
        value /= 1000;
    }
    return ( negv ? '-' : '' ) + Number( value ).toFixed( value < 10 && sfx != '' ? 1 : 0 ) + sfx;
};

// Implementation of Array.findIndex() and Array.find() functions for IE
if (Ext.isIE) {
    if (![].find) {
        Array.prototype.find = function(fn, scope) {
            var ret;
            if (!this.some(function(a, i) {
                    ret = a;
                    return fn.call(scope || this, a, i, this);
                })) return null;
            return ret;
        }
    }
    if (![].findIndex) {
        Array.prototype.findIndex = function(fn, scope) {
            var ret;
            if (!this.some(function(a, i) {
                    ret = i;
                    return fn.call(scope || this, a, i, this);
                })) return -1;
            return ret;
        }
    }
}

/* Cross-Browser Split 1.0.1
(c) Steven Levithan <stevenlevithan.com>; MIT License
An ECMA-compliant, uniform cross-browser split method */

var cbSplit;

// avoid running twice, which would break `cbSplit._nativeSplit`'s reference to the native `split`
if (!cbSplit) {

cbSplit = function (str, separator, limit) {
    // if `separator` is not a regex, use the native `split`
    if (Object.prototype.toString.call(separator) !== "[object RegExp]") {
        return cbSplit._nativeSplit.call(str, separator, limit);
    }

    var output = [],
        lastLastIndex = 0,
        flags = (separator.ignoreCase ? "i" : "") +
                (separator.multiline  ? "m" : "") +
                (separator.sticky     ? "y" : ""),
        separator = RegExp(separator.source, flags + "g"), // make `global` and avoid `lastIndex` issues by working with a copy
        separator2, match, lastIndex, lastLength;

    str = str + ""; // type conversion
    if (!cbSplit._compliantExecNpcg) {
        separator2 = RegExp("^" + separator.source + "$(?!\\s)", flags); // doesn't need /g or /y, but they don't hurt
    }

    /* behavior for `limit`: if it's...
    - `undefined`: no limit.
    - `NaN` or zero: return an empty array.
    - a positive number: use `Math.floor(limit)`.
    - a negative number: no limit.
    - other: type-convert, then use the above rules. */
    if (limit === undefined || +limit < 0) {
        limit = Infinity;
    } else {
        limit = Math.floor(+limit);
        if (!limit) {
            return [];
        }
    }

    while (match = separator.exec(str)) {
        lastIndex = match.index + match[0].length; // `separator.lastIndex` is not reliable cross-browser

        if (lastIndex > lastLastIndex) {
            output.push(str.slice(lastLastIndex, match.index));

            // fix browsers whose `exec` methods don't consistently return `undefined` for nonparticipating capturing groups
            if (!cbSplit._compliantExecNpcg && match.length > 1) {
                match[0].replace(separator2, function () {
                    for (var i = 1; i < arguments.length - 2; i++) {
                        if (arguments[i] === undefined) {
                            match[i] = undefined;
                        }
                    }
                });
            }

            if (match.length > 1 && match.index < str.length) {
                Array.prototype.push.apply(output, match.slice(1));
            }

            lastLength = match[0].length;
            lastLastIndex = lastIndex;

            if (output.length >= limit) {
                break;
            }
        }

        if (separator.lastIndex === match.index) {
            separator.lastIndex++; // avoid an infinite loop
        }
    }

    if (lastLastIndex === str.length) {
        if (lastLength || !separator.test("")) {
            output.push("");
        }
    } else {
        output.push(str.slice(lastLastIndex));
    }

    return output.length > limit ? output.slice(0, limit) : output;
};
cbSplit._compliantExecNpcg = /()??/.exec("")[1] === undefined; // NPCG: nonparticipating capturing group
cbSplit._nativeSplit = String.prototype.split;
} // end `if (!cbSplit)`

// for convenience... interferes with Ext.js interpretation of split, sorry
//String.prototype.split = function (separator, limit) {
//    return cbSplit(this, separator, limit);
//};
// End of Cross-Browser Split

String.prototype.trimToPix = function(length) {
    var tmp = trimmed = utils.sanitize(this);
    if (tmp.visualLength() > length)  {
        trimmed += "...";
        while (trimmed.visualLength() > length)  {
            tmp = tmp.substring(0, tmp.length-1);
            trimmed = tmp + "...";
        }
    }
    return trimmed;
};

String.prototype.trim = function () {
    return this.replace(/^\s*/, "").replace(/\s*$/, "");
};

String.prototype.visualLength = function() {
   var ruler = document.getElementById('string_ruler_unit');
   ruler.innerHTML = this;
   var ret = ruler.offsetWidth;
   ruler.innerHTML = '';
   return ret;
};

utils.ClearBrowserSelection = function() {
    try {
        var sel = window.getSelection ? window.getSelection() : document.selection;
        if (!sel) return;
        if (sel.removeAllRanges) sel.removeAllRanges();
        else 
            if (sel.empty) sel.empty(); 
    } catch(e) {}
};

utils.isNumeric = function(str) { return /^-?[0-9]+(\.[0-9]*)?[km]?$/i.test(str); };

utils.stringToNum = function(pos_str) {
    if (!pos_str) return;
    pos_str = pos_str.replace(/[, ]/g,'');
    if (pos_str.length < 1) return;
    var multiplier = 1;
    var last_char = pos_str.charAt(pos_str.length - 1).toUpperCase();
    if (last_char == 'K' || last_char == 'M') {
        pos_str = pos_str.substr(0, pos_str.length - 1);
        if (last_char == 'K') {
           multiplier = 1000;
        } else {
           multiplier = 1000000;
        }
    }
    var dec_part = 0;
    if (multiplier > 1) {
        var dec_pos = pos_str.indexOf('.');
        if (dec_pos != -1) {
            dec_part = Math.floor(parseFloat(pos_str.substr(dec_pos)) * multiplier);
            pos_str = pos_str.substr(0, dec_pos);
        }
    }
    return multiplier * parseInt(pos_str) + dec_part;
};

// Escape handles symbols innocuous from the point of view of HTML/HTTP
// but used internally in names, which are ususally provided by the user
// and as such are out of control. Encoding schema is according to SV-591
// and SV-1379:
// Prepend \ | , : [ ] with backslash, encode & ; # % as \hex-code, encode
// ' " = and space using standard %hex notation.
// Added a parameter 'symbs' allowing to process only a specific set of symbols
utils.escapeName = function(s, symbs) {
    var res = "";
    var parts = [];
    if (!symbs) parts = cbSplit(s, /([\]\[\\\|\'\"= ,:;&#%])/);
    else parts = cbSplit(s, symbs);
    for (var len = parts.length, i = 0; i < len;) {
        res += parts[i++];
        var sym = parts[i++];
        if (!sym) break;
        if (/[\]\[\\\|,:]/.test(sym))
            res += "\\" + sym;
        else {
            if (/['"= ]/.test(sym))
                res += "%";
            else
                res += "\\"; // ; & # %
            res += sym.charCodeAt(0).toString(16);
        }
    }
    return res;
};

utils.escapeTrackName = function(s, symbs) {
    var res = "";
    if (!symbs) symbs = /([\]\[\\\|,:=&;"#%+])/;
    var parts = cbSplit(s, symbs);
    for (var len = parts.length, i = 0; i < len;) {
        res += parts[i++];
        var sym = parts[i++];
        if (!sym) break;
        res += "\\";
        if (/[=&\\;"#%+]/.test(sym)) {
            res += sym.charCodeAt(0).toString(16);
        } else {
            res += sym;
        }
    }
    return res;
};

utils.unescapeName = function(s) {
    var parts = unescape(s).split("\\");
    var res = parts[0];
    for (var len = parts.length, i = 1; i < len; i++) {
        var part = parts[i];
        if (part.length == 0)
            res += "\\" + parts[++i];
        else if (/[\]\[\\\|,:]/.test(part.charAt(0)))
            res += part;
        else if (/^[0-9a-f]{2}.*/.test(part)) {
            res += String.fromCharCode(parseInt(part.slice(0,2), 16)) + part.slice(2);
        }
    }
    return res;
};


utils.sanitize = function(s) {
    return s.replace(/&/g, "&amp;").replace(/</g, "&lt;")
            .replace(/>/g, "&gt;").replace(/"/g, "&quot;");
};

// Draw a horizontal ruler with labels on top of it
// Parameters:
//     canvas  - canvas DOM element to draw in
//     base    - canvas base y coordinate, corresponds to ruler horizontal bar
//     start   - canvas x coordinate of ruler start
//     end     - canvas x coordinate of ruler end
//     mdStart - model start, i.e. start label number
//     mdEnd   - model end
utils.drawRuler = function(canvas, base, start, end, mdStart, mdEnd) {

    function logg() {
        //console.log.apply(null, arguments);
    }

    function getLabel(num) {
        return ""+num;
    }

    // scale - how one pixel is projected in model space
    function calcStep(ctx, scale, max_num) {
        function calcTickStep(base_step, scale) {
            var tick_step = base_step;
            [10, 5, 2].forEach(function(n) {
                var step = tick_step / n;
                //logg("n",n,"tick step ", step);
                if (step >= 1 && step / scale > 5) {
                    tick_step = step;
                    return;
                }
            });
            return tick_step;
        }

        var scale_sign = scale < 0 ? -1 : 1;
        scale = Math.abs(scale);

        var logOf10 = Math.log(10);

        var step = Math.pow(10, Math.ceil(Math.log(max_num)/logOf10));

        var maxLabelW = 2*ctx.measureText(""+max_num).width;
        logg("maxLabelW " + maxLabelW);
        var char_w = ctx.measureText("9").width;
        var comma_w = ctx.measureText(",").width;

        var base_step = step;
        if (step > 10.001) {
            var groups_n = 0;
            step = step * 10; // to compensate effect of the first iteration
            var max_label_w = maxLabelW;
            do
            {
                step = step / 10;
                var log = Math.ceil(Math.log(step+1)/logOf10);
                // number of comma-separated groups (111,222,333)
                groups_n = Math.ceil(log / 3);
                if (groups_n)  {
                // abbreviation correction coefficient for 100000 -> 100K etc.
                // should be coordinated with getLabel function
//                    var d_digits =  3 * groups_n - 2;
//                    max_label_w = maxLabelW - d_digits * char_w + (groups_n-1) * comma_w;
                    max_label_w = maxLabelW + (groups_n-1) * comma_w;
                    max_label_sym = scale * max_label_w;
                    logg("step "+step+" groups "+groups_n+" max_label_w " + max_label_w + " max_label_sym " + max_label_sym);
                }
            } while (groups_n  &&  step > max_label_sym*10);
            maxLabelW = max_label_w;
            base_step = step;
            // currently step has 10^X form, lets check if we can choose
            // a smaller step in a form  K * 10^(X-1), where K = 2, 5
            // this adjusment does not affect labels size
            if (step > max_label_sym * 5)    {
                base_step = step / 10;  // 10^(X-1)
                step = step / 5;        // 2 * 10^(X-1)
            }
            else if (step > max_label_sym * 2)   {
                base_step = step / 10;  // 10^(X-1)
                step = step / 2;        // 5 * 10^(X-1)
            }
        }
        var tick_step = calcTickStep(base_step, scale);
        logg("step "+step+" base_step "+base_step+" tick_step "+tick_step);
        logg("step "+step/scale+"px base_step "+base_step/scale+"px");
        return [step, base_step, tick_step, scale_sign];
    }

    var pixLen = end - start + 1;
    var mdLen = mdEnd - mdStart + 1;

    // Here we adjust the pixel region so that the ticks are in the middle
    // of a displayed letter
    var pix_per_letter = Math.abs(pixLen/mdLen);
    if (pix_per_letter > 1) {
        start = Math.round(start + pix_per_letter/2);
        end   = Math.round(end - pix_per_letter/2);
        pixLen = end - start + 1;
    }

    var scale = (mdLen-1) / pixLen;
    logg("scale "+scale);
    var ctx = canvas.getContext("2d");
    ctx.font = "11px Helvetica";
    var res = calcStep(ctx, scale, Math.max(mdStart, mdEnd));
    var mdStep = Math.max(1, res[0]) * res[3],
        mdBaseStep = Math.max(1, res[1]) * res[3],
        mdTickStep = Math.max(1, res[2]) * res[3];
    var labelFreq = mdStep/mdTickStep,
        baseFreq  = mdBaseStep/mdTickStep;
    ctx.textAlign = "center";
    ctx.lineWidth = 1;
    ctx.beginPath();
    var mdStartOffset = Math.ceil(mdStart%mdTickStep);
    if (mdStartOffset) mdStartOffset -= mdTickStep;
    var md_labelOffset = Math.ceil(mdStart%mdStep);
    if (md_labelOffset) md_labelOffset -= mdStep;
    md_labelOffset -= mdStartOffset;
    logg("mdStartOffset", mdStartOffset,
         "md_labelOffset", md_labelOffset);
    var tickStart = Math.ceil(start-mdStartOffset/scale),
        tickEnd = end + 0.5,
        tickStride = mdTickStep/scale,
        tickHeight = 1,
        tickBaseHeight = 2,
        tickLabelHeight = 3,
        tickBase = base + 0.5,
        ticks = pixLen / tickStride + 1,
        labelBase = tickBase - tickLabelHeight-4;

    // draw bar
    ctx.beginPath();
    ctx.moveTo(start+0.5, tickBase);
    ctx.lineTo(end+0.5, tickBase);
    ctx.stroke();
    // draw end ticks
    ctx.moveTo(start+0.5, tickBase - tickLabelHeight);
    ctx.lineTo(start+0.5, tickBase + tickLabelHeight);
    ctx.stroke();
    ctx.moveTo(end+0.5, tickBase - tickLabelHeight);
    ctx.lineTo(end+0.5, tickBase + tickLabelHeight);
    ctx.stroke();

    // draw end labels
    ctx.textAlign = "left";
    var labelText = getLabel(mdStart);
    ctx.fillText(labelText, start, labelBase);
    var px_leftLimit = start + ctx.measureText(labelText).width;
    ctx.textAlign = "right";
    var labelText = getLabel(mdEnd);
    ctx.fillText(labelText, end, labelBase);
    var px_rightLimit = end - ctx.measureText(labelText).width;

    // draw inner ticks and labels
    ctx.textAlign = "center";
    var mdStride = (mdEnd - mdStart) / (ticks-1);
    for (var px_pos = tickStart, tick = 0; px_pos < tickEnd; px_pos += tickStride, ++tick) {
        var x = Math.ceil(px_pos) + 0.5;
        var isLabel = (mdTickStep*tick+md_labelOffset) % mdStep == 0;
        var isBase = (mdTickStep*tick+md_labelOffset) % mdBaseStep == 0;
        var tick_up = isLabel ? tickLabelHeight : 0;
        var tick_down = isLabel ? tickLabelHeight : (isBase ? tickBaseHeight : tickHeight);
        ctx.moveTo(x, tickBase - tick_up);
        ctx.lineTo(x, tickBase + tick_down);
        ctx.stroke();
        if (isLabel) {
            var labelText = getLabel(mdStart-mdStartOffset + Math.round(tick*mdStride));
            var width = ctx.measureText(labelText).width;
            if (x - width/2 > px_leftLimit && x + width/2 < px_rightLimit)
                ctx.fillText(labelText, x, labelBase);
        }
    }
}

utils.clearAria = function(div) {  // clear aria-owns
    Ext.query('[aria-owns]:not([aria-owns=""])', div).forEach(function(el) {
        el.removeAttribute('aria-owns');
    });
    Ext.query('[aria-readonly="false"]', div).forEach(function(el) {
        el.removeAttribute('aria-readonly');
    });
}

/**
*  MD5 (Message-Digest Algorithm)
*  http://www.webtoolkit.info/
**/

utils.MD5 = function (string) {
    function RotateLeft(lValue, iShiftBits) {
        return (lValue<<iShiftBits) | (lValue>>>(32-iShiftBits));
    }
    function AddUnsigned(lX,lY) {
        var lX4,lY4,lX8,lY8,lResult;
        lX8 = (lX & 0x80000000);
        lY8 = (lY & 0x80000000);
        lX4 = (lX & 0x40000000);
        lY4 = (lY & 0x40000000);
        lResult = (lX & 0x3FFFFFFF)+(lY & 0x3FFFFFFF);
        if (lX4 & lY4) {
            return (lResult ^ 0x80000000 ^ lX8 ^ lY8);
        }
        if (lX4 | lY4) {
            if (lResult & 0x40000000) {
                return (lResult ^ 0xC0000000 ^ lX8 ^ lY8);
            } else {
                return (lResult ^ 0x40000000 ^ lX8 ^ lY8);
            }
        } else {
            return (lResult ^ lX8 ^ lY8);
        }
    }
    function F(x,y,z) { return (x & y) | ((~x) & z); }
    function G(x,y,z) { return (x & z) | (y & (~z)); }
    function H(x,y,z) { return (x ^ y ^ z); }
    function I(x,y,z) { return (y ^ (x | (~z))); }
    function FF(a,b,c,d,x,s,ac) {
        a = AddUnsigned(a, AddUnsigned(AddUnsigned(F(b, c, d), x), ac));
        return AddUnsigned(RotateLeft(a, s), b);
    };
    function GG(a,b,c,d,x,s,ac) {
        a = AddUnsigned(a, AddUnsigned(AddUnsigned(G(b, c, d), x), ac));
        return AddUnsigned(RotateLeft(a, s), b);
    };
    function HH(a,b,c,d,x,s,ac) {
        a = AddUnsigned(a, AddUnsigned(AddUnsigned(H(b, c, d), x), ac));
        return AddUnsigned(RotateLeft(a, s), b);
    };
    function II(a,b,c,d,x,s,ac) {
        a = AddUnsigned(a, AddUnsigned(AddUnsigned(I(b, c, d), x), ac));
        return AddUnsigned(RotateLeft(a, s), b);
    };
    function ConvertToWordArray(string) {
        var lWordCount;
        var lMessageLength = string.length;
        var lNumberOfWords_temp1=lMessageLength + 8;
        var lNumberOfWords_temp2=(lNumberOfWords_temp1-(lNumberOfWords_temp1 % 64))/64;
        var lNumberOfWords = (lNumberOfWords_temp2+1)*16;
        var lWordArray=Array(lNumberOfWords-1);
        var lBytePosition = 0;
        var lByteCount = 0;
        while ( lByteCount < lMessageLength ) {
            lWordCount = (lByteCount-(lByteCount % 4))/4;
            lBytePosition = (lByteCount % 4)*8;
            lWordArray[lWordCount] = (lWordArray[lWordCount] | (string.charCodeAt(lByteCount)<<lBytePosition));
            lByteCount++;
        }
        lWordCount = (lByteCount-(lByteCount % 4))/4;
        lBytePosition = (lByteCount % 4)*8;
        lWordArray[lWordCount] = lWordArray[lWordCount] | (0x80<<lBytePosition);
        lWordArray[lNumberOfWords-2] = lMessageLength<<3;
        lWordArray[lNumberOfWords-1] = lMessageLength>>>29;
        return lWordArray;
    };
    function WordToHex(lValue) {
        var WordToHexValue="",WordToHexValue_temp="",lByte,lCount;
        for (lCount = 0;lCount<=3;lCount++) {
            lByte = (lValue>>>(lCount*8)) & 255;
            WordToHexValue_temp = "0" + lByte.toString(16);
            WordToHexValue = WordToHexValue + WordToHexValue_temp.substr(WordToHexValue_temp.length-2,2);
        }
        return WordToHexValue;
    };
    function Utf8Encode(string) {
        string = string.replace(/\r\n/g,"\n");
        var utftext = "";
        for (var n = 0; n < string.length; n++) {
            var c = string.charCodeAt(n);
            if (c < 128) {
                utftext += String.fromCharCode(c);
            }
            else if((c > 127) && (c < 2048)) {
                utftext += String.fromCharCode((c >> 6) | 192);
                utftext += String.fromCharCode((c & 63) | 128);
            }
            else {
                utftext += String.fromCharCode((c >> 12) | 224);
                utftext += String.fromCharCode(((c >> 6) & 63) | 128);
                utftext += String.fromCharCode((c & 63) | 128);
            }
        }
        return utftext;
    };
    var x=[];
    var k,AA,BB,CC,DD,a,b,c,d;
    var S11=7, S12=12, S13=17, S14=22;
    var S21=5, S22=9 , S23=14, S24=20;
    var S31=4, S32=11, S33=16, S34=23;
    var S41=6, S42=10, S43=15, S44=21;

    string = Utf8Encode(string);
    x = ConvertToWordArray(string);
    a = 0x67452301; b = 0xEFCDAB89; c = 0x98BADCFE; d = 0x10325476;

    for (k=0;k<x.length;k+=16) {
        AA=a; BB=b; CC=c; DD=d;
        a=FF(a,b,c,d,x[k+0], S11,0xD76AA478);
        d=FF(d,a,b,c,x[k+1], S12,0xE8C7B756);
        c=FF(c,d,a,b,x[k+2], S13,0x242070DB);
        b=FF(b,c,d,a,x[k+3], S14,0xC1BDCEEE);
        a=FF(a,b,c,d,x[k+4], S11,0xF57C0FAF);
        d=FF(d,a,b,c,x[k+5], S12,0x4787C62A);
        c=FF(c,d,a,b,x[k+6], S13,0xA8304613);
        b=FF(b,c,d,a,x[k+7], S14,0xFD469501);
        a=FF(a,b,c,d,x[k+8], S11,0x698098D8);
        d=FF(d,a,b,c,x[k+9], S12,0x8B44F7AF);
        c=FF(c,d,a,b,x[k+10],S13,0xFFFF5BB1);
        b=FF(b,c,d,a,x[k+11],S14,0x895CD7BE);
        a=FF(a,b,c,d,x[k+12],S11,0x6B901122);
        d=FF(d,a,b,c,x[k+13],S12,0xFD987193);
        c=FF(c,d,a,b,x[k+14],S13,0xA679438E);
        b=FF(b,c,d,a,x[k+15],S14,0x49B40821);
        a=GG(a,b,c,d,x[k+1], S21,0xF61E2562);
        d=GG(d,a,b,c,x[k+6], S22,0xC040B340);
        c=GG(c,d,a,b,x[k+11],S23,0x265E5A51);
        b=GG(b,c,d,a,x[k+0], S24,0xE9B6C7AA);
        a=GG(a,b,c,d,x[k+5], S21,0xD62F105D);
        d=GG(d,a,b,c,x[k+10],S22,0x2441453);
        c=GG(c,d,a,b,x[k+15],S23,0xD8A1E681);
        b=GG(b,c,d,a,x[k+4], S24,0xE7D3FBC8);
        a=GG(a,b,c,d,x[k+9], S21,0x21E1CDE6);
        d=GG(d,a,b,c,x[k+14],S22,0xC33707D6);
        c=GG(c,d,a,b,x[k+3], S23,0xF4D50D87);
        b=GG(b,c,d,a,x[k+8], S24,0x455A14ED);
        a=GG(a,b,c,d,x[k+13],S21,0xA9E3E905);
        d=GG(d,a,b,c,x[k+2], S22,0xFCEFA3F8);
        c=GG(c,d,a,b,x[k+7], S23,0x676F02D9);
        b=GG(b,c,d,a,x[k+12],S24,0x8D2A4C8A);
        a=HH(a,b,c,d,x[k+5], S31,0xFFFA3942);
        d=HH(d,a,b,c,x[k+8], S32,0x8771F681);
        c=HH(c,d,a,b,x[k+11],S33,0x6D9D6122);
        b=HH(b,c,d,a,x[k+14],S34,0xFDE5380C);
        a=HH(a,b,c,d,x[k+1], S31,0xA4BEEA44);
        d=HH(d,a,b,c,x[k+4], S32,0x4BDECFA9);
        c=HH(c,d,a,b,x[k+7], S33,0xF6BB4B60);
        b=HH(b,c,d,a,x[k+10],S34,0xBEBFBC70);
        a=HH(a,b,c,d,x[k+13],S31,0x289B7EC6);
        d=HH(d,a,b,c,x[k+0], S32,0xEAA127FA);
        c=HH(c,d,a,b,x[k+3], S33,0xD4EF3085);
        b=HH(b,c,d,a,x[k+6], S34,0x4881D05);
        a=HH(a,b,c,d,x[k+9], S31,0xD9D4D039);
        d=HH(d,a,b,c,x[k+12],S32,0xE6DB99E5);
        c=HH(c,d,a,b,x[k+15],S33,0x1FA27CF8);
        b=HH(b,c,d,a,x[k+2], S34,0xC4AC5665);
        a=II(a,b,c,d,x[k+0], S41,0xF4292244);
        d=II(d,a,b,c,x[k+7], S42,0x432AFF97);
        c=II(c,d,a,b,x[k+14],S43,0xAB9423A7);
        b=II(b,c,d,a,x[k+5], S44,0xFC93A039);
        a=II(a,b,c,d,x[k+12],S41,0x655B59C3);
        d=II(d,a,b,c,x[k+3], S42,0x8F0CCC92);
        c=II(c,d,a,b,x[k+10],S43,0xFFEFF47D);
        b=II(b,c,d,a,x[k+1], S44,0x85845DD1);
        a=II(a,b,c,d,x[k+8], S41,0x6FA87E4F);
        d=II(d,a,b,c,x[k+15],S42,0xFE2CE6E0);
        c=II(c,d,a,b,x[k+6], S43,0xA3014314);
        b=II(b,c,d,a,x[k+13],S44,0x4E0811A1);
        a=II(a,b,c,d,x[k+4], S41,0xF7537E82);
        d=II(d,a,b,c,x[k+11],S42,0xBD3AF235);
        c=II(c,d,a,b,x[k+2], S43,0x2AD7D2BB);
        b=II(b,c,d,a,x[k+9], S44,0xEB86D391);
        a=AddUnsigned(a,AA);
        b=AddUnsigned(b,BB);
        c=AddUnsigned(c,CC);
        d=AddUnsigned(d,DD);
    }
    var temp = WordToHex(a)+WordToHex(b)+WordToHex(c)+WordToHex(d);

    return temp.toLowerCase();
}


} //if (window.NCBIGBUtils)
/*  $Id: tooltip.js 46857 2021-11-09 20:58:12Z borodine $
 * File Description: Ext.ToolTip extention
*/
if (!window.NCBIGBObject || !NCBIGBObject.ToolTip) {
Ext.tip.QuickTipManager.init(true);// , {showDelay: 1000, mouseOffset: [-15, 15]});

if (Ext.supports.Touch && Ext.versions.core.version.indexOf('6.0.') == 0) {
    Ext.define('Ext.tip.ToolTip', {
        override: 'Ext.tip.ToolTip',
        setTarget: function(target) {
            var me = this,
                t = Ext.get(target),
                tg;
            if (me.target) {
                tg = Ext.get(me.target);
                if (Ext.supports.Touch) me.mun(tg, 'tap', me.onTargetOver, me);
                me.mun(tg, {
                    mouseover: me.onTargetOver,
                    mouseout: me.onTargetOut,
                    mousemove: me.onMouseMove,
                    scope: me
                });
            }
            me.target = t;
            if (t) {
                if (Ext.supports.Touch)  me.mon(t, { tap: me.onTargetOver, scope: me });
                me.mon(t, {
                    mouseover: me.onTargetOver,
                    mouseout: me.onTargetOut,
                    mousemove: me.onMouseMove,
                    scope: me
                });
            }
            if (me.anchor) {
                me.anchorTarget = me.target;
            }
        }
    });
}

Ext.define('NCBIGBObject.ToolTip', {
    extend: 'Ext.tip.ToolTip',
    hideDelay: 300,
    showDelay: 1000,
    draggable: true,
    pinned: false,
    autoScroll: false,
    adjustWidth: false,
    pinCallbackFn: Ext.emptyFn,
    header: {
        baseCls: 'x-panel-header',
        padding: 0,
        titlePosition: 99
    },
    listeners: {
        beforehide: function() {
            return !(this.insideTT || this.pinned);
        }
    },
    initComponent: function() {
        if (this.pinnable) delete this.pinnable;
        this.callParent(arguments);
        if (this.isPinnable())
            this.addTool([{type: 'unpin', scope: this, hidden: false, callback: this.togglePin}]);
     },
/*
     afterRender: function() {
         this.callParent(arguments);
//         this.body.setStyle('position', 'absolute');
    },*/
    isPinnable: function() { return this.pinnable != false; },

    doWidthAdjusting: function() {
        if (!this.adjustWidth || !this.isVisible()) return; 
        var maxW = 490;
            minW = this.getWidth();
        if (maxW < minW || maxW - minW < 50 || minW <= this.minWidth) return;

        var delta = 0,
            curH = this.getHeight(),
            minH = this.setWidth(maxW).getHeight();

        while (minH < curH) {
            delta = maxW - minW;
            delta = delta >> 1;
            curH = this.setWidth(minW + delta).getHeight();
            if (minH < curH++) minW += delta;
            else maxW -= delta;

            if (delta < 3) break;
        }
        this.setWidth(minW + delta);
    },

    onMouseOver: function(e) {
        e.stopEvent();
        this.insideTT = true;
    },

    onMouseOut: function(e) {
        e.stopEvent();
        this.insideTT = e.within(this.el, true, true);
        if (!(this.insideTT || this.pinned)) this.hide();
    },

    onShow: function() {
        this.callParent(arguments);
        this.getEl().on('mouseover', this.onMouseOver, this);
        this.getEl().on('mouseout', this.onMouseOut, this);
        this.doWidthAdjusting();
    },

    togglePin: function(o, tool) {
        if (!this.pinned) {
            this.pinCallbackFn();
            tool.setType('pin');
            if (!this.isVisible()) this.show();
        } else {
            tool.setType('unpin');
        }
        this.autoHide = !(this.pinned = !this.pinned);
    },
    update: function(html) {
        if (this.adjustWidth) delete this.width; // = 'auto';
        NCBIGBObject.ToolTip.superclass.update.call(this, html);
        this.updateLayout();
        this.doWidthAdjusting();
        if (this.adjustHeight || (!this.isHidden() && this.getWidth() >= 500)) {
            this.setScrollY(25);
            var scrl = this.getScrollY();
            if (scrl && scrl < 25) this.setHeight(this.getHeight() + scrl + 5);
            else this.setScrollY(0);
        }
    }
});
}

/*  $Id: base.js 46757 2021-09-28 22:06:28Z borodine $
 *
 * File Description: base classes 
 *
 */


if (!window.NCBIGBClass) {
Ext.define('NCBIGBClass.View', {
    constructor: function(type, app) {
        this.clear();
        this.m_Type = type;
        this.m_App = app;
        this.m_ReqNum = 0;
        NCBIGBUtils.m_NextViewIdx = NCBIGBUtils.m_NextViewIdx || 0;
        this.m_Idx = NCBIGBUtils.m_NextViewIdx++;
        this.m_View = null;
    },

    isAlignment: function() { return false; },
    isLoading: function() { return this.m_Loading; },
    isGraphic: function() { return false; },
    isPanorama: function() { return false; },
    getBodyDiv: function() { return this.bodyDiv || (this.bodyDiv = Ext.get(this.m_DivId)); },
    getBodyWidth: function() { return this.getBodyDiv().getWidth(); },
    getHeight: function() { return this.m_Height; },
    getMostRightPix: function() { return 0; },
    getSpacerHeight: function() { return 4; },
    getWidth: function() { return this.m_Width; },
    getXY: function() { return this.m_View.getEl().getXY(); },
    moveTo: function(vis_from, vis_len) {},
    ping: function(a) { this.m_App.ping(a); },
    pingClick: function(a, e) { this.m_App.pingClick(a, e); },
    refresh: function() {},
    updateTitle: function() {},
    updateTracks: function() {},
    updateStatusInfo: function() {},

    clear: function() {
        this.m_Loading =  false;
        this.m_ScrollPix = this.m_Width = this.m_Height = 0;
        this.m_Theme = this.m_DivId = this.m_FromCgi = this.m_View = this.m_Color = this.m_Spacer = this.m_Locator = null;
        delete this.m_lastParams;
    },

    createCloseTBar: function() {
        return  {
            id:'close', qtip:'Close View', scope:this,
            handler:function(e, target, panel) {
                this.m_App.viewIsClosing(this);
                if(panel.view) { panel.view.remove(); }
            }
        };
    },

    destroy: function() {
        if (this.m_Locator) {
            this.m_Locator.remove();
        }
        if(this.m_View) {
            if (this.m_Spacer) {
                this.m_View.ownerCt.remove(this.m_Spacer, true);
            }
            this.m_View.ownerCt.remove(this.m_View, true);
        }
        this.clear();
    },

    remove: function() {
        this.m_App.removeView(this);
        this.destroy();
    }
});


/* locator */
NCBIGBClass.Locator = function(view, color, resizable) {
    this.m_View = view; // view
    var panBody = view.m_App.m_Panorama.getBodyDiv();

    var tpl_pan = new Ext.Template('<div class="pan-bar" style="background-color:#{color};position:absolute;top:0px;left:0px;"></div>');
    var tpl = new Ext.Template('<div style="position:absolute;top:17px;left:0px;width:1px;z-index:5;"><div class="locator_rect"></div></div>');
    this.scroller = tpl_pan.insertFirst(panBody, {color:color}, true);
    this.rectangle = tpl.insertFirst(panBody, null, true);
    this.rectangle.setHeight(this.getPanoramaHeight() - 2);

    this.scroller.on({
        'pointerdown'  : this.onPointerDown,
        'contextmenu':  this.onContextMenu,
        scope: this
    });

    if (resizable) {
        (new Ext.Template('<div class="left-resizer"></div>')).insertAfter(this.rectangle.getLastChild());
        this.rectangle.getLastChild().on({
            'pointerdown': this.onPointerDown,
            scope: this
        });
        (new Ext.Template('<div class="right-resizer"></div>')).insertAfter(this.rectangle.getLastChild());
        this.rectangle.getLastChild().on({
            'pointerdown': this.onPointerDown,
            scope: this
        });
        this.setColor({value: color});
    }

};

NCBIGBClass.Locator.prototype = {
    getPanoramaWidth: function() { return this.m_View.m_App.getPanoramaWidth(); },
    getPanoramaHeight: function() { return this.m_View.m_App.getPanoramaHeight(); },
    getRight:  function(local) { return this.rectangle.getRight(local); },
    getWidth:  function(contentWidth) { return this.rectangle.getWidth(contentWidth); },
    setWidth: function(width) { this.rectangle.setWidth(width); },
    setHeight:  function(h) { this.rectangle.setHeight(h - 17); },
    getLeft:  function(local) { return this.rectangle.getLeft(local); },
    setLeft:  function(pos) {
        var panorama_width = this.getPanoramaWidth();
        var bar_pos = pos;
        if(pos + 24 > panorama_width)
            bar_pos = panorama_width - 24;

        this.scroller.setLeft(bar_pos);
        this.rectangle.setLeft(pos);
    },

    setColor:  function(cpicker) {
        this.m_View.m_Color = cpicker.value;
        this.scroller.setStyle('background-color', '#' + cpicker.value);
        this.m_View.m_View.body.setStyle('borderLeft', '#' + cpicker.value + ' 1px solid');
    },


    remove:  function() {
        this.scroller.remove();
        this.rectangle.remove();
    },

    onPointerDown: function(e) {
        if (this.m_ContextMenu) this.m_ContextMenu.destroy();
        if ((e.pointerType == 'mouse' && e.button) || this.m_View.m_App.m_Panorama.m_Loading) return;
        this.m_View.m_App.m_Panorama.m_Locator = this;
        this.downX = Math.round(e.getX());
        this.offsetX = this.getLeft(true) - this.downX;
        this.m_Action = '1-0-R';
        switch (e.getTarget().className) {
            case 'left-resizer':
                this.moveHandler = function(newX) {
                    var new_left = Math.min(this.getRight(true) - 2, Math.max(0, newX));
                    this.setWidth(this.getRight(true) - new_left);
                    this.setLeft(new_left);
                }
                break;
            case 'right-resizer':
                this.offsetX = this.getRight(true) - this.downX;
                this.moveHandler = function(newX) {
                    var new_right = Math.max(this.getLeft(true) + 2,  Math.min(newX, this.getPanoramaWidth()));
                    this.setWidth(new_right - this.getLeft(true));
                }
                break;
            case 'pan-bar':
                if (e.pointerType == 'touch') this.m_deferredContext = Ext.defer(this.onContextMenu, 2000, this);
                this.moveHandler = function(newX) {
                    var new_left = Math.min(this.getPanoramaWidth() - this.getWidth(),  Math.max(0, newX));
                    this.setLeft(new_left);
                    this.m_Action = '1-0-D';
                }
                break;
            default: 
                delete this.offsetX;
                return;
        }
        if (Ext.isIE) e.stopPropagation(); else e.stopEvent();
        this.m_View.m_App.eventHandlers('on', this);
    },

    onPointerUp: function(e) {
        if (this.m_deferredContext) this.m_deferredContext = clearTimeout(this.m_deferredContext);
        e.stopPropagation();
        this.m_View.m_App.eventHandlers('un', this);
        if (typeof this.offsetX != 'undefined') {
            this.m_View.syncToLocator(this.m_Action.slice(-1) == 'D' ? this.m_View.m_VisLenSeq : 0);
            this.m_View.pingClick(this.m_Action);
            delete this.offsetX;
        }
    },

    onPointerMove: function(e) {
        if (typeof this.offsetX == 'undefined') return;
        e.stopPropagation();
        //e.stopEvent();
        NCBIGBUtils.ClearBrowserSelection();
        var majicSensivity = 3;
        if (e.pointerType == 'touch' && this.m_deferredContext) {
            if (Math.abs(this.downX - e.getX()) > majicSensivity) {
                clearTimeout(this.m_deferredContext);
                this.m_deferredContext = 0;
            } else return;
        }
        this.moveHandler(Math.round(e.getX()) + this.offsetX);
    },

    onContextMenu: function(e) {
        var menu = new Ext.menu.Menu();
        delete this.m_Action;
        if (!this.m_deferredContext) {
            e.stopEvent();
        } else {
            this.m_deferredContext = 0;
            this.m_ContextMenu = menu;
        }
        menu.add([{
            text: 'Bring to front', scope: this, handler:function() {
                this.rectangle.insertBefore(Ext.get('pan-holder' + this.m_View.m_App.m_Idx));
                this.scroller.insertAfter(this.rectangle);
            }}, {
            text: 'Sent to back', scope: this, handler:function() {
                this.m_View.m_App.m_Panorama.bodyDiv.insertFirst(this.rectangle);
                this.scroller.insertAfter(this.rectangle);
            }}, {
            text: 'View color change',  menu: new Ext.menu.ColorPicker({listeners: { 'select': this.setColor, scope: this }})}
        ]);
        if (this.m_View.m_App.m_Panorama.modifyLocatorCtxMenu) this.m_View.m_App.m_Panorama.modifyLocatorCtxMenu(menu); 
        menu.showAt(this.scroller.getLeft(false), this.scroller.getBottom(false));
    }
};

Ext.define('NCBIGBClass.Panorama', {
    extend: 'NCBIGBClass.View',
    m_PrevXY: null,
    m_vLine: {setX: Ext.emptyFn, setStyle: Ext.emptyFn},

    constructor: function(app) {
        this.callParent(['panorama', app]);
        this.m_TopOffset = 18; //sm_HolderSize
        this.m_DivId = 'panorama_' + app.m_DivId + app.m_Idx;
        this.m_selectionDivId = 'sel_' + this.m_DivId;
        this.m_View = this.m_App.addView({ 
            html:'<div id="' + this.m_DivId + '" class="panorama_div"><div id="' + this.m_selectionDivId
            + '" class="panoramaSelection" style="z-index:1;">'
            + '</div><div id="pan-holder' + this.m_App.m_Idx + '" class="pan-holder"/>'
        });
        Ext.get(this.m_DivId).on({ scope: this,
            pointerdown: this.onPointerDown,
            pointerup:   this.onPointerUp,
            pointermove: this.onPointerMove,
            contextmenu: this.onContextMenu,
            pointerleave: this.onPointerLeave
        });

        //if  embedded=panorama then do not allow selection
        if (!this.m_App.m_AllViewParams.match("embedded=panorama")) 
            this.m_PanoramaSelection = new NCBIGBClass.PanoramaSelection(this.m_DivId,this.m_selectionDivId, [0,0]);
    },

    onContextMenu: Ext.emptyFunc,

    onPointerLeave:  function(e) {
        this.m_vLine.setStyle('display', 'none');
    },

    isPanorama: function() { return true; },
    
    toPix: function(x) { return x * this.m_Width / this.m_App.getDataLength(); },// sequence 2 screen (panorama)
    toSeq: function(x) { return Math.round(this.m_App.getDataLength() * x / this.m_Width); },
    seq2Pix: function(seq_pos) { return this.toPix(seq_pos); },
    seq2PixScrolled: function(seq_pos) { return this.seq2Pix(seq_pos) + this.m_ScrollPix; },// For panorama m_ScrollPix is always zero, but anyway
    pix2Seq: function(pix_pos) { return this.toSeq(pix_pos); },// screen 2 sequence (gview)
    getMostRightPix: function() { return this.getWidth(); },
    
    onPointerDown: function(e) {
        if (this.m_ContextMenu) this.m_ContextMenu.destroy();
        if (e.pointerType == 'mouse' && e.button) return
        var tID = e.target.id;
        if (this.m_PanoramaSelection && tID.search('scroller') < 0 && tID.search('resizer') < 0){
            this.m_ResizeAction = true;
            this.m_XY = e.getXY();
            if (e.pointerType != 'mouse') {
                this.m_deferredContext = Ext.defer(this.onContextMenu, 2000, this, null);
            } else this.m_deferredContext = 0;
            this.m_XFinal_Selection = this.m_XY[0] + 1;
            this.m_PanoramaSelection.resize(this.m_XY[0], this.m_XFinal_Selection);
            var el =  Ext.get(this.m_selectionDivId);
           if (el) el.applyStyles("display:inline;");
           e.stopEvent();
        }
    },

    onPointerMove: function(e) {
        var magicSensivity = 5; 
        var x = e.getX();
        if (!this.m_ResizeAction || !this.m_PanoramaSelection) {
            this.m_vLine.setStyle('display', 'block');
            this.m_vLine.setX(x);
            return;
        }
        this.m_vLine.setStyle('display', 'none');

        if (this.m_deferredContext) {
            if (Math.abs(x - this.m_XY[0]) <= magicSensivity) return;   
            clearTimeout(this.m_deferredContext);
            this.m_deferredContext = 0;
        }
        this.m_XFinal_Selection = x;
        this.m_PanoramaSelection.resize(this.m_XY[0], x);
        this.m_Selection = true;
        e.stopEvent();
    },

    onPointerUp: function(e) {
        var el =  Ext.get(this.m_selectionDivId);
        if (!el || !el.isVisible() || !this.m_PanoramaSelection) return;
        if (this.m_deferredContext) {   
            clearTimeout(this.m_deferredContext);
            this.m_deferredContext = 0;
        }
        this.m_ResizeAction = false;
        if (Math.abs(this.m_XFinal_Selection - this.m_XY[0]) < 5) {
            el.applyStyles("display:none;");
            return;
        } else {
//            if(this.m_App.m_Views[1]){
                this.m_Locator.setLeft(el.getLeft() - Ext.get(this.m_DivId).getLeft());
                this.m_Locator.setWidth(el.getWidth());
                this.m_Locator.m_View.syncToLocator();
                el.applyStyles("width:0px;display:none;");
                this.pingClick('1-0-1');
//            }
        }
        e.stopPropagation();
    },

    onClick: function(e) {
        var el =  Ext.get(this.m_selectionDivId);
        el.setWidth(0);
        el.hide();
        this.m_ResizeAction = false;
        
    },


    onResize: function(e) {},
    
    createMarkerElem: function(marker) {
        var elem = Ext.get(this.m_DivId);
        var create_params = marker.getCreateParams(true,this.m_Idx);
        var marker_elem = create_params.template.append(elem, create_params.options, true);
        marker_elem.setTop(this.m_TopOffset + 3);
        return marker_elem;
    }
    

});

//////////////////////////////////////////////////////////////////////////
NCBIGBClass.PanoramaSelection = function(panoramaDivId, selectionDivId, range) {
    this.m_selectionDivId = selectionDivId;
    this.m_panoramaDivId = panoramaDivId;
    this.m_Resizing = false;
    this.range = range;
    this.element = Ext.get(this.m_selectionDivId);
};

NCBIGBClass.PanoramaSelection.prototype = {

    pageToViewX: function(x) {
        var div_xy = Ext.get(this.m_panoramaDivId).getXY();
        return x - div_xy[0];
    },

    resize: function(xSelection1, xSelection2) {
        var x1 = this.pageToViewX(xSelection1);
        var x2 = this.pageToViewX(xSelection2);
        var el = this.element;
        el.setLeft(x1 <= x2 ? x1 : x2);
        el.setWidth(Math.abs(x2 - x1));
    }
  
};
}




/*  $Id: sv_base.js 45927 2021-01-13 17:59:54Z borodine $
 *
 * Authors:  
 *
 * File Description: Sequence Viewer modification of base classes
 *
 */

Ext.define('SeqView.View', {
    extend: 'NCBIGBClass.View',

    createChooseColorBtn: function() {
        var colorBarTPL = new Ext.Template('<div id="view_color_{idx}" style="width:15px;height:15px;background-color:#{color};border:1px black solid;">&nbsp;</div>');
        if (!this.m_Color || this.m_Color.length == 0) {
            var cp = new Ext.ColorPalette();
            this.m_Color = cp.colors[Math.floor( Math.random() * cp.colors.length)];
        }
        return {
                text: colorBarTPL.apply({idx:this.m_Idx, color: this.m_Color}),
                tooltip:'View Color',
                scope: this,
                handler: function() {this.pingClick('1-2-0');},
                menu:new Ext.menu.ColorPicker({listeners: {'select': function(cm, color) {
                    Ext.get('view_color_'+this.m_Idx).setStyle('background-color', '#'+color);
                    this.m_Color = color;
                    if (this.m_Locator) { this.m_Locator.setColor(cm); }
                },scope:this}})
        };
    },

    createCloseTBar: function() {
        return  {
            id:'close', qtip:'Close View', scope:this,
            handler:function(e, target, panel) {
                this.m_App.viewIsClosing(this);
                var view = panel.view;
                if(view) { view.remove(); }
            }
        };
    },

    getHeight: function(noruler) { return this.m_Height - (!noruler ? 0 : ((this.m_bottomRulerHeight || 0) + (this.m_topRulerHeight || 0))); },
    setGeneMode: function(){}

});

SeqView.Locator = NCBIGBClass.Locator;

Ext.define('SeqView.Panorama', {
    extend: 'NCBIGBClass.Panorama',


    updateView: function(callback) {
        this.m_Width = this.getBodyWidth();
        if (this.m_Width <= 0) {
            // the view is hidden, don't do anything now.
            return;
        }
        this.callbackOnLoad = callback || Ext.emptyFn;
        var url = this.m_App.m_CGIs.Panorama;
        var params = {theme:'NCBI Overview', id: this.m_App.GI, assm_context: this.m_App.m_AssmContext, width: this.m_Width, affinity_token: this.m_App.defaultAffinityToken};
        Ext.apply(params, this.m_App.m_GraphicExtraParams);
        this.m_App.addURLParams(params);
        if (Ext.Object.equals(this.m_lastParams || {}, params)) {
            this.callbackOnLoad(); 
            return; //no need to re-render
        }
        this.m_lastParams = params;
        this.m_Loading = true;
        this.m_App.AjaxRequest({url: url, context: this, data: params,
             success:this.checkJobStatus, error:this.loadFailure});
        delete this.m_lastParams.appname;
    }, 

    checkJobStatus: function(data, text, res) {
        var from_cgi = SeqView.decode(data);
        if (from_cgi.job_status) {
            if (from_cgi.job_status == 'failed') {
                this.loadFailure(null, from_cgi.error_message);
            } else if(from_cgi.job_status == 'canceled') {
                this.loadFailure(null, 'Job canceled');
            } else {
                var url = this.m_App.m_CGIs.Panorama + '?job_key=' + from_cgi.job_id
                Ext.defer(this.m_App.AjaxRequest, 2000,this,[{url: url, context: this,
                        success: this.checkJobStatus, error: this.loadFailure}]);
            }
        } else {
            if (from_cgi.success === false) {
                this.loadFailure(null, from_cgi.msg);
            } else {
                // If the img_url begins with ? it contains only parameters for ncfetch, so prepend ncfetch URL
                // This is a way to provide reliable URL resolution for embedding. SV-1760
                if (from_cgi.img_url && from_cgi.img_url.charAt(0) == '?') {
                    from_cgi.img_url = this.m_App.m_CGIs.NetCache + from_cgi.img_url;
                }
                this.m_FromCgi = from_cgi;
                this.m_Height = this.m_FromCgi.img_height + this.m_TopOffset + 1;
    
                var the_div = Ext.get(this.m_DivId);
                var img_el = the_div.first('img');
                if (img_el) {
                    var d = Ext.getDom(img_el);
                    d.src = from_cgi.img_url;
                } else {
                    var tpl = new Ext.Template('<img aria-label = "panorama" src="{img_url}">'),
                        img_el = tpl.append(the_div, from_cgi,true);

                }
                the_div.setHeight(this.m_Height);
                this.m_View.updateLayout();
                this.m_vLine = the_div.appendChild(new Ext.Element(document.createElement('div')));
                this.m_vLine.set({style: 'position:absolute;border:1px solid lightgray;border-right-width:0; top:0; width:0px; height:'
                    + this.m_Height + 'px; display: none;'});

                this.m_App.updateMarkersSize(this); // update markers
                this.m_App.updateMarkersPos(this); // update markers
                this.m_Loading = false; // loaded
                this.m_App.notifyViewLoaded(this);
                this.callbackOnLoad();
                if( !this.m_App.m_DialogShown ){
                    this.m_App.resizeIFrame();
                }
            }
        }
    },

    loadFailure: function(data, text, res) {
        Ext.MessageBox.show({title:'Image loading error', msg: text, buttons: Ext.MessageBox.OK, icon: Ext.MessageBox.INFO});
    },

    refresh: function() {
        var app = this.m_App;
        this.updateView(function() { app.forEachView(function(v) {app.updateLocator(v);}); });
    },

    onContextMenu: function(e) {
        var menu = new Ext.menu.Menu();
        if (this.m_deferredContext) {
            clearTimeout(this.m_deferredContext);
            this.m_deferredContext = 0;
            this.m_ContextMenu = menu;
        } else {
            this.m_XY = e.getXY();
            e.preventDefault();  // this prevents the default contextmenu to open in Firefox (linux)
            e.stopPropagation();
        }

        this.pingClick('1-1');
        var x_pos = this.m_XY[0] - Ext.get(this.m_DivId).getX()
        menu.add({text:'Set New Marker At Position', iconCls:'xsv-markers', scope:this,  
                handler:function() { this.m_App.newMarkerDlg(this, x_pos); this.pingClick('1-1-1');} 
        });
        menu.add('-');
        if (this.m_App.m_Origin) {
            menu.add({text:'Reset Sequence Origin', iconCls:'xsv-origin', scope:this, 
                handler:function() { this.m_App.clearOrigin(); this.pingClick('1-1-2'); }
            });
        } else {
            menu.add({text:'Set Sequence Origin At Position', iconCls:'xsv-origin', scope:this, 
                handler:function() { this.m_App.setOrigin(this, x_pos); this.pingClick('1-1-2'); }
            });
        }
        var nvMarker = this.m_App.getMarkersInfo().findMarkerByName(SeqView.MarkerNav);
        if (nvMarker) {
            menu.add({text:'Remove Navigation Marker', iconCls: 'xsv-marker-remove', scope:this,  
                handler:function() { nvMarker.deleteMarker(); this.pingClick('1-1-3'); }
            });
        }
        if (menu.items.length>0) { menu.showAt(this.m_XY); }
    }
});
/*  $Id: selection.js 47423 2023-03-28 16:11:16Z borodine $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Authors:  Vlad Lebedev, Maxim Didenko, Victor Joukov
 *
 * File Description: Mouse-over selection
 *
 */

Ext.define('SeqView.SelectionToolTip', {
    extend: 'NCBIGBObject.ToolTip',
    dismissDelay: 0,
    minWidth: 200,
    autoScroll: true,
    cls: 'SeqViewerApp',

    initComponent : function(){
        this.gview = this.selection.m_View;
        this.area = this.selection.area;
        this.pinnable = (this.area.type & SeqView.AreaFlags.NoPin) == 0;
        if (this.area.type & SeqView.AreaFlags.Track) { //track ToolTip
            this.pinnable = false;
            this.track = this.gview.getTrackByArea(this.area);
            if (this.track.hub_name && this.area.tooltip.indexOf(this.track.hub_name + '</a>') < 0)
                this.area.tooltip += ' <br>Track Hub: <a href='
                    + this.track.hub_url + ' target=\"_blank\">' + this.track.hub_name + '</a>';
            if (this.track.data_provider.indexOf('_rmt' > 0)) {
                var url = this.track.annots.find(function(a) { return a.search(/(https?|ftp):\/\//) == 0; });
                if (url && this.area.tooltip.indexOf('Track source:') < 0)
                    this.area.tooltip += ' <br>Track source: <a href='
                    + url + ' target=\"_blank\">' + url.split('/').pop() + '</a><br>';
            }
        }
        this.target = this.selection.element;

        this.callParent(arguments);

        if (this.track) {
            this.setMinWidth(280);
            this.setMaxWidth(400);
            this.addTool([
/*                { type: 'gear',  callback: this.trackTTHandler, tooltip: 'Modify Track Settings', stat: '8-1', scope: this,
                  hidden: !(this.track.check_boxes || this.track.choice_list)},*/
                { type: 'restore', callback: this.trackTTHandler, 'data-qtip': 'Track standalone view', stat: '8-4', scope: this }
//                { type: 'help', callback: this.trackTTHandler, qtip: 'Track legend', stat: '8-2', scope: this },
//                { type: 'close', hidden: this.gview.m_App.m_NoConfDlg, callback: this.trackTTHandler, tooltip: 'Hide (deactivate) Track', stat: '8-3', scope: this }
            ]);
            this.addDocked(
                { xtype:'toolbar', dock: 'bottom', border: false,
                    items: [
                    {xtype: 'button', text: 'Track Legend',
                     hidden: !this.track.legend_text, icon: SeqView.base_url + 'images/info.png', scope: this, stat: '8-2',
                     handler: this.buttonsHandler},
                    {xtype: 'button', text: 'View in CGV', iconCls: 'xsv-cgv',
                     hidden: !this.area.links || this.area.links[0].label != 'CGV', scope: this, stat: '8-cgv',
                     handler: this.buttonsHandler},
                    '->',
                    {xtype: 'button', text: 'Configure', qtip: 'Modify Track Settings', tooltip: 'Modify Track Settings', style: {pointerEvents: 'all'},
                     hidden: !this.track.hasOpts, iconCls: 'xsv-config', scope: this, stat: '8-1', handler: this.buttonsHandler},
                    {xtype: 'button', text: 'Download', tooltip: 'Download Track Data', icon: SeqView.base_url + 'images/downloadtrk.png',
                     scope: this, stat: '8-dwnld', hidden: this.gview.m_App.m_noTrackData || !SeqView.TM.isDownloadable(this.track),
                     handler: this.buttonsHandler}
                    ]}
             );
        } else {
            this.title = this.area.title || this.selection.title;
            if (this.isPinnable()) {
                this.pinCallbackFn = function() {
                    this.down('#search').show();
                    this.gview.pingClick('6-1');
                }
                this.addTool([
                    { type: 'search', callback: this.showFeat, itemId: 'search', hidden: true, scope: this },
                    { type: 'magnify', callback: this.loadAreaRange, scope: this },
                    { type: 'up', itemId: 'up', callback: this.collapseHandler, hidden: true, scope: this },
                    { type: 'down', itemId: 'down', callback: this.collapseHandler, hidden: true, scope: this }]);
                }
        }
    },

    buttonsHandler: function(btn, e) {
        this.gview.pingClick(btn.stat);
        switch (btn.stat){
            case '8-1': SeqView.TM.modifyTrackDetails(this.gview, this.track); break;
            case '8-2': SeqView.showHelpDlg('legends/#' + this.track.legend_text); break;
            case '8-cgv': window.open(this.area.links[0].link + this.gview.m_VisFromSeq + '-' + (this.gview.m_VisFromSeq + this.gview.m_VisLenSeq));break;
            default: this.gview.downloadTblFile([this.track]);
        }
        this.destroy();
    },

    collapseHandler: function(o, tool, e) {
        this.collapsedTT = tool.type == 'up';
        this.setTTContent();
        this.show(this.getXY());
        this.gview.pingClick('6-4');
    },

    trackTTHandler: function(o, tool, e) {
        var app = this.selection.m_View.m_App;
        this.hide();
        this.gview.pingClick(tool.stat);
        switch (tool.type) {
            case 'close':
                this.track.shown = false;
                if (this.gview.m_SelectedSig && !this.gview.m_featMarkers)
                    SeqView.TM.Common.updateSeqViewApp(SeqView.TM.processTracksInfo(app.m_Config.TrackConfig), app);
                else { 
                    this.gview.hideTrack(this.track.cfgIdx);
                    this.destroy();
                }
                break;
            case 'gear': SeqView.TM.modifyTrackDetails(this.selection.m_View, this.track); break;
//            case 'help': SeqView.showHelpDlg('legends/#' + this.track.legend_text); break;
            case 'restore':
                var str = '';
                Ext.each(app.m_Config.TrackConfig, function(trk) {
                    if (trk.shown && trk.key == "sequence_track") {
                        str += SeqView.TM.tracksArrayToString(trk);
                        return false;
                    }
                });
                this.gview.m_graphTracks.some(function(t){
                    if (this.track.cfgIdx == t.cfgIdx) {
                        app.getLinkToThisPageURL('full', str + t.tParams);
                        return false;
                    }
                }, this);
                break;
        }
    },

    onShow: function() {
        if (this.gview.m_ContextMenu && this.gview.m_ContextMenu.isVisible()) {
            this.gview.m_ContextMenu.hide();
        }
        this.callParent(arguments);
    },

    setTTContent: function(tt) {
        var collapsedToolTips = 'NCBI/SV/ToolTips/collapse';
        if (typeof this.collapsedTT == 'undefined' && typeof localStorage !== 'undefined')
            this.collapsedTT = !(localStorage[collapsedToolTips] === undefined);
        else
            if (typeof localStorage !== undefined) {
                if (!this.collapsedTT) localStorage.removeItem(collapsedToolTips);
                else localStorage.setItem(collapsedToolTips, '');
            }

        tt = tt || this.area.tooltip;
        if (!tt) {
            this.mask('Loading...');
//            this.update('<div style="white-space:nowrap;">Data loading...</div>');
            this.adjustWidth = false;
            return;
        }
        if (this.isMasked() && this.getEl()) this.unmask();

        if (typeof tt == 'string') {
            this.update(tt);
            return;
        }

        if (tt.title && !this.track) {this.setTitle(''); this.setTitle(tt.title)};

        var content = tt.mainInfo + (tt.downloadLinks ? '<br/>' + tt.downloadLinks : '');

        if (tt.basicLinks) content += '<br/><b>Links & Tools</b><br/>' + tt.basicLinks;

        var tUp, tDown;
        if (tt.extraLinks) {
            tUp = this.header.down('#up');
            tDown = this.header.down('#down');
            if (tUp) tUp.hide();
            if (tDown) tDown.hide();
            if (!this.collapsedTT) {
               if (!tt.basicLinks) content += '<br/><b>Links & Tools</b>';
               content += '<br/>' + tt.extraLinks;
               tDown = false;
            } else tUp = false;
        }
        content += tt.regionBtns ? tt.regionBtns.html.replace(/selectionID/g, this.id) : '';

        this.update(content);
        if (tUp) tUp.show();
        if (tDown) tDown.show();
    },

    onRender: function() {
        if (!this.area.tooltip) {
            SeqView.App.simpleAjaxRequest(this.selection.ajaxCfg);
        }
        this.callParent(arguments);
    },

    onHide: function() {
        this.callParent(arguments);
        if (this.pinned) {
            delete this.selection.area.tooltip.selectionTT;

        }
        this.destroy();
    },


    afterRender : function() {
        this.callParent(arguments);
        this.setTTContent()
    },

    showFeat: function(tt) {
        var area = this.gview.m_FromCgi.areas.find(function(a) { return a.signature == tt.area.signature; })
        if (area) {
            this.area = area;
            var view_start = -this.gview.m_ScrollPix;
            if (this.area.bounds.r - 10 < view_start) {
                this.gview.scrollViewTo(-(this.area.bounds.l-10));
            } else  if (this.area.bounds.l+10 > this.gview.getBodyWidth() + view_start) {
                var new_pos = view_start + (this.area.bounds.r - (view_start + this.gview.getBodyWidth()));
                this.gview.scrollViewTo(-(new_pos+10));
            }
            this.highlightFeat();
        } else {
            this.loadAreaRange(false);
        }
        this.gview.pingClick('6-2');
    },

    highlightFeat: function() {
        this.gview.blinkAreas({blinkData: {area: this.gview.updateArea(this.area)}});
    },
    
    loadAreaRange: function(magnify) {
        var from = this.area.range[0];
        var len = this.area.range[1] - this.area.range[0] + 1;
        if (magnify) {
            var la = Math.floor(len/10);
            from = Math.max(0, from - la);
            len += 2 * la;
            this.gview.pingClick('6-3');
        } else {
            if (len < this.gview.m_VisLenSeq) {
                from = Math.max(0, from - (this.gview.m_VisLenSeq - len)/2);
                len = this.gview.m_VisLenSeq;
           }
        }
        if (this.gview.m_VisFromSeq != from || this.gview.m_VisLenSeq != len) {
            this.gview.startImageLoading(from, len, {from_ui: true, doMove: true, cbScope: this, callback: this.highlightFeat});
        } else {
            this.highlightFeat();
        }
    },

    highlightToolTip: function() {
        this.getEl().addCls('xsv-tooltip_highlight');
    },
    unHighlightToolTip : function(){
        this.getEl().removeCls('xsv-tooltip_highlight');
    }

});


SeqView.createSelectionElement = function(view, parent_elem, area, tpl) {
    var bounds = area.bounds;
    var element = tpl.append(parent_elem, {idx: Ext.id()}, true);

    var left  = Math.max(0, view.getFlip() ? bounds.r : bounds.l);
    var right  = Math.min(view.m_Width, view.getFlip() ? bounds.l : bounds.r);
    var the_width = Math.min(view.m_Width, right - left) + 2;
    var the_top = bounds.t - 1;
    var the_left  = left - 1;// + view.m_ScrollPix;
/*
    if (area.type && (area.type & SeqView.AreaFlags.Track)) {
        the_left = 0;
        the_width = view.bodyDiv.getWidth();
    }*/
    if (area.tname) {
        the_top = bounds.t;
        the_width = bounds.r - bounds.l;
    }
    element.setLeft(the_left); // the output from cgi seems a bit off. That's why this offsets. Or is it browser-specific? :)
    element.setTop(the_top);
    element.setHeight(bounds.b-bounds.t);
    element.setWidth(the_width);
    return element;
}


SeqView.Selection = function(view, parent_div, areas, event, offsetX) {
    var app = view.m_App,
        flags = SeqView.AreaFlags,
        signatures = [],
        disable_caching = true;
    this.m_View = view;
    this.delayTimeout = null;
    this.saved_event = null;
    this.areas = areas;
    this.area = areas.length ? areas[areas.length - 1] : areas;

    var trk = (this.area.type & flags.Track) ? this.m_View.getTrackByArea(this.area) : null;
//    if (event.ctrlKey && trk) console.log(trk.cfgIdx, trk);

    for (var i = 0; i != areas.length; i++) {
        signatures.push( areas[i].signature );
        if (areas[i].feats_count) {
            this.area.merged_feats = areas[i].feats_count;
        }
        disable_caching &= Boolean(areas[i].type & flags.NoCaching);
    }
    if (disable_caching && areas.length) this.area.type |= flags.NoCaching;

    var tpl = new Ext.Template('<div id="selection_id_{idx}" class="over_' +
            ((this.area.type & (flags.NoHighlight | flags.Track | flags.Sequence)) !== 0 ? 'no' : '') +
            'selection sv-drag sv-dblclick"/>');

    this.element = SeqView.createSelectionElement(view, parent_div, this.area, tpl);

    this.element.on({ scope: this,
        'click':       (!trk || trk.hasOpts) ? this.clickTest : Ext.emptyFn,
        'contextmenu': this.onContextMenu,
        'dblclick':    this.interOnDblClick,
        'mouseout':    this.remove
    });

    var area = this.area,
        params = {appname: app.m_AppName, bases_per_pix: 1/view.m_bpPix, color: app.m_Config.Options.curr_color, id: app.GI, objinfo: 1, signatures: signatures.join(',')};
    if (typeof this.area.tooltip == 'string' && this.area.type == '0') {
        this.area.textTT = this.area.tooltip;
        delete this.area.tooltip;
    }
    if (!area.tooltip) {
        if (area.type & (flags.Sequence | flags.Ruler | flags.TooltipEmbedded)) area.tooltip = "Not provided";
        if (trk && !area.tooltip && !area.signature.match(";")) area.tooltip = trk.help;
        if (area.tooltip && area.tooltip.match("Not provided")) return;
    }
    if (!area.title) area.title = area.range ? (area.range[0] + 1) + ".." + (area.range[1] + 1) : area.name;

    if ((area.type & flags.Track) && area.signature.match(";")){
        var gijunk = area.signature.split(";");
        params.signatures = gijunk.pop().split("|")[1];
    }

    if (!(area.type & flags.NoCaching)) {
        area.tooltip = area.tooltip || this.m_View.m_App.findSelectionToolTip(params.signatures);
    } else if (!(area.type & flags.Track) && !(area.type & flags.TooltipEmbedded)) {
        delete area.tooltip;
    }
    
    var tt = area.tooltip;
    if (tt) { // use from cache
        if (tt.selectionTT && tt.selectionTT.pinned) {
           this.m_ToolTip = tt.selectionTT;
           this.m_ToolTip.selection = this;
           this.m_ToolTip.setTarget(this.element);
           this.m_ToolTip.highlightToolTip();
        }
        else this.m_ToolTip = new SeqView.SelectionToolTip({ selection: this });
    } else {
        var trk = this.m_View.getTrackByArea(area);
        var str = '[id:' + trk.id + ',' + 'key:' + trk.key;
        Ext.each(['name', 'data_key', 'display_name', 'dbname'],
            function(p) { if (trk[p]) str += ',' + p + ':' + NCBIGBUtils.escapeTrackName(trk[p]); });
        Ext.each(trk.annots, function(annot, idx){
            str += (idx == 0 ? ',annots:' : '|') + NCBIGBUtils.escapeTrackName(annot);
        });
        params.track = str + ']';
        if (trk.affinity_token) params.affinity_token = trk.affinity_token;

        var key = (app.m_Key || '') + ',' + (trk.seq_data_key || '')
        if (key != ',') params.key = key.replace(/^,|,$/,'');
        if (app.m_DepthLimit) params.depthlimit = app.m_DepthLimit;
        if (app.m_BamPath) params.bam_path = app.m_BamPath;
        if (app.m_SRZ) params.srz = app.m_SRZ;
        if (app.m_AssmContext) params.assm_context = app.m_AssmContext;
        if (area.merged_feats) { params.merged_feats = area.merged_feats; }
        if (offsetX !== null && (area.type & flags.NoCaching)) {
            this.element.on('mousemove',function(e) { 
                this.ajaxCfg.data.seq_pos = view.getSeqPos(e.getX() - offsetX);
            }, this);
            params.seq_pos = view.getSeqPos(event.getX() - offsetX);
            if (app.getFlip()) params.flip = 1;
        }
        this.ajaxCfg = { url: trk.pURL, data: params, context: this,
                        success: this.checkJobStatus };
        this.m_ToolTip = new SeqView.SelectionToolTip({ selection: this });
    }
    // 'Tap' event
    if (event && event.type == 'pointerup') {
        this.m_ToolTip.showAt(event.getXY());
        this.m_ToolTip.setY(event.getY()); //Workaround to fix ExtJS glitch with Y position when a page is scrolled
    }
};


SeqView.Selection.prototype = {
    movePix: function(delta) {
        if (this.removed) return;
        var new_left = this.element.getLeft(true) - delta;
        this.element.setLeft(new_left);
    },

    remove: function() {
        if (this.removed) return;
        if (this.m_ToolTip) {
            if (this.m_ToolTip.pinned) this.m_ToolTip.unHighlightToolTip();
            if (this.area.tooltip) this.area.tooltip.selectionTT = this.m_ToolTip;
            this.m_ToolTip.setTarget();
        }
        this.element.remove();
        this.removed = true;
    },

    onClick: function(e) {
        if (SeqView.scrolled) {
            clearTimeout(SeqView.scrolled);
            delete SeqView.scrolled;
            return;
        }
        if (this.area.type & SeqView.AreaFlags.Track)
            SeqView.TM.modifyTrackDetails(this.m_View, this.m_View.getTrackByArea(this.area));
        else
            this.m_View.changeSelectedSig(this.area, e.ctrlKey);
        if (!this.m_ToolTip.pinned) this.m_ToolTip.destroy();
    },

    onDblClick: function(e) {
        //for ruler and track names area boundaries are ignored and zoom ins are based on coordinates
        if (this.area.type & (SeqView.AreaFlags.Track | SeqView.AreaFlags.Ruler)) return;
        if (this.area.type & (SeqView.AreaFlags.ZoomInOnDblClck)) {
            this.m_View.onDblClick(e);
            return;
        }
        var range = this.area.range;
        var new_from  = (range[0] < range[1]) ? range[0] : range[1];
        var new_to    = (range[0] < range[1]) ? range[1] : range[0];
        var new_len = new_to - new_from + 1; // before offset
        var l_offset = new_len / 100 * 6;  // calc. 5% offset
        new_from -= l_offset;
        new_to += l_offset;
        new_len = new_to - new_from + 1; // after offset
        this.m_View.startImageLoading(new_from, new_len, {from_ui: true});
    },

    fireClick: function(e) {
        clearTimeout(this.delayTimeout);
        this.delayTimeout = null;
        if (!this.DblClickHit)
            this.onClick(this.saved_event);
        this.DblClickHit = false;
        this.saved_event = null;
    },

    clickTest: function(e) {
        if (this.delayTimeout) {
            if (e.type == 'click') {
                this.DblClickHit = false;
                clearTimeout(this.delayTimeout);
                this.delayTimeout = null;
                this.saved_event.type = 'dblclick';
                this.saved_event.browserEvent = 'Event dblclick';
                this.onDblClick(this.saved_event);
                this.saved_event = null;
            }
        } else {
            if (e.type == 'click' ){
                this.saved_event = new Ext.EventObjectImpl(e);
                var that = this;
                this.delayTimeout = setTimeout(function() {that.fireClick(e);}, 400);
            }
        }
    },

    interOnDblClick: function(e) {
        this.DblClickHit = true;
        this.onDblClick(e);
    },

    onContextMenu: function(e) {
        e.stopEvent();
        if (this.m_ToolTip) {
            this.m_ToolTip.destroy();
            this.m_ToolTip == null;
        }
        var range = this.area.range;
        if (range === undefined) return;
        var seq_pos = this.area.strand == 1 ? range[0]: range[1];
        var that = this;
        var areas;
        if (this.area.signature.indexOf('fake|') == -1) areas = this.areas;
        this.m_View.showContextMenu(e.getXY(), seq_pos/*,
            function(menu, download_menu) {
                if (areas) that.addFeatureLinks(menu, seq_pos, areas[0]);
                that.addDownloadLinks(download_menu, areas);
            }*/
        );
    },

    addFeatureLinks: function(menu, seq_pos, area) {
        var submenuVT = new Ext.menu.Menu();
        var menuVT = menu.add('-', {text: 'Views & Tools', iconCls: 'x-tbar-loading', disabled: true, menu: submenuVT})[1];
        var view = this.m_View;

        var addSubmenu = function(links) {
            menuVT.setIconCls('xsv-views_tools')
            if (links.length) menuVT.enable();

            Ext.each(links, function(l) {
                submenuVT.add({text: l.label, url: l.link,
                    handler: function(conf) { view.pingClick('2-1-VTlink'); window.open(conf.url); }});
            });
        }
        if (area.links) {
            addSubmenu(area.links);
            return;
        }
        var params = this.createCGIParams(area.signature);
        params.link = 1;
        this.m_View.m_App.AjaxRequest({
            url: this.m_View.m_App.m_CGIs.Graphic, context:this, data: params,
            success: function(data, txt, res) {
                var links = SeqView.decode(data);
                for (var i = 0; i < links.length; i++) {
                    var url = links[i].link;
                    for (var j = 0; j < i; j++) {
                        if (links[j].link == url) {
                            links.splice(i--, 1);
                            break;
                        }
                    }
                }
                if (this.m_View.m_App.m_Embedded !== false)
                    for (var i = 0; i < links.length; i++) links[i].link = links[i].link.replace(/^\//, SeqView.webNCBI);
                area.links = links;
                addSubmenu(links);
            }
        }); // Ext.Ajax

    },

    addDownloadLinks: function(menu, areas) {
        if (!menu) return;
        var params = this.createCGIParams();
        var val = "";
        for (var i = 0; i < areas.length; i++) {
            if (val) val += ",";
            val += areas[i].signature;
        }
        params.ids = val;
        params.link_type = "download";
        params.link = 1;
        var view = this.m_View;
        view.m_App.AjaxRequest({
            url: view.m_App.m_CGIs.Graphic,
            context: this, data: params,
            success: function(data, txt, res) {
                var links = SeqView.decode(data);
                for (var i = 0; i < links.length; i++) {
                    menu.add({text:links[i].label, scope:this, link: links[i].link,
                            handler: function(self) {
                                var fn = self.link.split('&filename=')[1].split('&', 1)[0];
                                view.pingClick('2-1-9-featFA');
                                view.m_App.requestFile(view, self.link.split('.cgi',2)[1], '&filename=' + fn);
                            }
                    });
                }
            }
        }); // Ext.Ajax
    },

    createCGIParams: function(signature) {
        var params = {};
        if (signature)
            params.id = signature;
        if (this.m_View.m_App.m_ViewerContext)
            params.viewer_context = this.m_View.m_App.m_ViewerContext;
        if (this.m_View.m_App.m_Key)
            params.key = this.m_View.m_App.m_Key;
        if (this.m_View.m_App.m_DepthLimit)
            params.depthlimit = this.m_View.m_App.m_DepthLimit;
        if (this.m_View.m_App.m_SRZ)
            params.srz = this.m_View.m_App.m_SRZ;
        if (this.m_View.m_App.m_BamPath)
            params.bam_path = this.m_View.m_App.m_BamPath;
        return params;
    },

    checkJobStatus: function(data, txt, req) {
        var objs = SeqView.decode(data);
        if (objs.job_status && (objs.job_status != 'failed' || objs.job_status != 'canceled')) {
            this.ajaxCfg.data = {job_key: objs.job_id};
            Ext.defer(SeqView.App.simpleAjaxRequest, 1000, this, [this.ajaxCfg]);
            return;
        }
        var obj = objs[0];
        if (obj.error && this.area.textTT) {
            obj.text = this.area.textTT;
            obj.title = this.area.title;
            delete obj.error;
        }
        if (!obj) obj = {text: this.m_View.getTrackByArea(this.area).help};
        var links = obj.links || [];
        var tip = obj.text || obj.error;
        if (objs.length > 1) { // multiple features
            tip += '<br>' + (objs[1].text || objs[1].error);
            links = links.concat(objs[1].links || []);
        }
        var app = this.m_View.m_App;
        if (app.m_preprocessorTT) {
            var prs = obj.object_id.split('-');
            app.m_preprocessorTT(links,
                {acc: app.m_Config.SeqInfo.id, title: obj.title,// GI: app.m_Config.SeqInfo.gi,
                start_pos: parseInt(prs[1], 16) + 1, end_pos: parseInt(prs[2], 16) + 1});
        }
        var pmids = tip.split('PMID:');
        tip = pmids[0];
        for (var i = 1; i < pmids.length; i++) {
            var pmid = parseInt(pmids[i]);
            if (!isNaN(pmid)) {
                tip += '<a href="https://www.ncbi.nlm.nih.gov/pubmed/' + pmid + '" target="_blank">PMID:';
                pmids[i] = pmids[i].replace(pmid, pmid + '</a>');
            }
            tip += pmids[i];
        }

        if (obj.exts) this.area.exts = obj.exts;
        if (obj.title) this.area.title = obj.title;

        var tt = {mainInfo: tip, basicLinks: '', extraLinks: '', downloadLinks: '', title: this.area.title};

        var linksSet = new Set();
        for (var j = 0, jlen = links.length; j < jlen; j++) {
            var lnk = links[j];
            var path = (lnk.type == 'Download' ? SeqView.base_url : SeqView.webNCBI).slice(0, -1);
            var str = '<td align="right" valign="top" nowrap><b>';
            var link_name = (lnk.name && lnk.name.length > 0) ? lnk.name + ': ' : '';
            var link_key = window.JSON.stringify(lnk);;
            if (linksSet.has(link_key))
                continue;
            linksSet.add(link_key);
            str += link_name.replace(/ /g, '&nbsp;');
            str += '</b></td><td valign="top">';
            var anchor = lnk.html || '';
            if (!anchor) {
                for (var k = 0; k < lnk.links.length; k++) {
                    var linrec = lnk.links[k];
                    var found = false;
                    for (var l = 0; l < k; l++) {
                        if (lnk.links[l].link == linrec.link){
                            found = true;
                            break;
                        }
                    }
                    if (found) continue;

                    if (k > 0) anchor += '<br>';
                    var label = linrec.label.replace(/ /g, '&nbsp;');
                    anchor += '<a href="'
                    if (linrec.link.search(/^https?:/i) && linrec.link.search(/^ftp:/i))
                        anchor += path + (linrec.link.charAt(0) == '/' ? '' : '/') ;
                    anchor += linrec.link + '" target="_blank">' + label + '</a>';
                }
            }
            if(lnk.type == 'Download') anchor = anchor.replace(/sequence.cgi\?/g, 'sequence.cgi?netcache=0&');
            str += anchor + '</td>';
            if(anchor) {
                switch (lnk.type) {
                    case 'Basic': tt.basicLinks += "<tr>" + str + "</tr>"; break;
                    case 'Extra': tt.extraLinks += "<tr>" + str + "</tr>"; break;
                    case 'Download': tt.downloadLinks += "<tr>" + str + "</tr>"; break;
                }
            }
        }
        var ttID = this.m_ToolTip.id,
            rData = [];
            rStr = '';
        Ext.each((obj.unaligned_regions || []).concat(obj.insertions || []), function (unr) {
            var idx = rData.push(unr) - 1,
                l = unr.length - 1,
                frwd = unr.orientation == 'forward';
            if (unr.to == undefined) unr.to = unr.from + unr.length - 1;
            unr.range = (1 + (frwd ? unr.from : unr.to)) + '..' + (1 + (!frwd ? unr.from : unr.to));
            id = idx + '_seq_id_' + ttID;
            rStr += '<tr><td align="right" valign="top" nowrap><b>' + unescape(unr.label || 'Insertion') + '</b> (' + unr.range + '):&nbsp;'
                 + '</b></td><td valign="top">' + (unr.to - unr.from + 1) + (unr.polya ? '&nbsp;(polyA)' : '') + '</td><td>&nbsp;&nbsp;&nbsp;';
            if (!unr.seq_id && unr.length < 10) rStr += unr.sequence;
            else
                rStr +=  '<button id=\"' + id + '\" title="View ' + (unr.seq_id ? 'unaligned fragment' : 'insertion') + '"'
                     + 'onclick="SeqView.clickButtonTT(\'selectionID\',' + idx + ')" '
                     + 'style=" font-size: 11px; border-radius: 4px; border-width: 1px; background-color: whitesmoke;">'
                     + 'View</button>';
            rStr += '</td></tr>';
        });
        if (rStr) {
            tt.regionBtns = {html: '<br/><b>Unaligned Regions&nbsp;</b><table border="0" cellpadding="0" cellspacing="0"><tbody>'
                + rStr + '</tbody></table>', unr: rData};
            this.adjustHeight = true;
        }
        if (tt.basicLinks) tt.basicLinks = '<table border="0" cellpadding="0" cellspacing="0"><tbody>' + tt.basicLinks + '</tbody></table>';
        if (tt.extraLinks) tt.extraLinks = '<table border="0" cellpadding="0" cellspacing="0"><tbody>' + tt.extraLinks + '</tbody></table>';
        if (tt.downloadLinks) tt.downloadLinks = '<table border="0" cellpadding="0" cellspacing="0"><tbody>' + tt.downloadLinks + '</tbody></table>';
        this.area.tooltip = tt;
        if (!(this.area.type & SeqView.AreaFlags.NoCaching))
            this.m_View.m_App.saveSelectionToolTip(this.ajaxCfg.data.signatures, tt);

        if (this.m_ToolTip.isVisible()) this.m_ToolTip.setTTContent(tt);

        this.m_View.m_App.ping({"jsevent": "mouseover","sv-event":"Object_ToolTip"});
    },

    showPrime: function(unr, params) {
        var color = 'cornflowerblue';
        if (unr.sequence) {
            SeqView.TM.showSeqRegion('Insertion (' + unr.range + ')', '<mark style="background-color:' + color + '">' + unr.sequence + '</mark>');
            return;
        }
        var urStart = (unr.label.search(/5|N/) == 0),
            from = unr.from,
            to = unr.to,
            reverse = unr.orientation == 'reverse',
            length = Math.abs(to - from),
            tail = Math.min(unr.length - (to - from) - 1, 500);
        if (length > 2000) {
            if (urStart) from = to - 2000;
            else to = from + 2000;
        }
        var area = this.area,
            trk = this.m_View.getTrackByArea(area);
        var str = '[id:' + trk.id + ',' + 'key:' + trk.key;
        Ext.each(['name', 'data_key', 'display_name', 'dbname'],
            function(p) { if (trk[p]) str += ',' + p + ':' + NCBIGBUtils.escapeTrackName(trk[p]); });
        Ext.each(trk.annots, function(annot, idx){
            str += (idx == 0 ? ',annots:' : '|') + NCBIGBUtils.escapeTrackName(annot);
        });
        params = params || {url: trk.pURL, context: this,
            data: {objinfo: 1, id: unr.seq_id, unalignedregion: true, from: from, to: to,
                alignedtail: tail, reverse: reverse, tailtype: urStart ? 'start' : 'end', track: str + ']'}};
        if (this.m_View.m_App.m_AssmContext) params.data.assm_context = this.m_View.m_App.m_AssmContext;
        if (this.m_View.m_App.getFlip()) params.data.flip = true;
        if (trk.pURL.indexOf('_rmt.cgi') > 0) params.data.signatures = area.signature;
        if (trk.affinity_token) params.data.affinity_token = trk.affinity_token;
        if (this.m_View.m_Flip != params.data.reverse) urStart = ! urStart;
        params.success = function(data){
            switch (data.job_status) {
                case 'submitted': case 'running': case 'pending':
                    params.data.job_id = data.job_id;
                    Ext.defer(this.showPrime, 1000, this, [unr, params]);
                    break;
                default:
                    var up = data.unaligned_part,
                        ap = data.aligned_part;
                    SeqView.TM.showSeqRegion('Unaligned Region (' + unr.orientation + ' strand)',
                        '<p style="word-wrap:break-word; margin: 3px;">' + (urStart ? '' : ap.sequence) + '<mark style="background-color:' + color + '">' + up.sequence + '</mark>' + (urStart ? ap.sequence : '') + '</p>');
            }
        }
        params.error = function(d, t){
            Ext.MessageBox.alert('Error', d);
        }
        this.m_View.m_App.AjaxRequest(params);
    }
};

SeqView.clickButtonTT = function(id, idx) {
    var selTT = Ext.getCmp(id);
    if (selTT) selTT.selection.showPrime(selTT.area.tooltip.regionBtns.unr[idx]);
};



Ext.define('SeqView.SelectedRangeToolTip', {
    extend: 'NCBIGBObject.ToolTip',
    anchor: 'top',
    anchorToTarget: true,
    initComponent: function() {
        this.gview = this.selection.m_View;
        this.target = this.selection.m_tt_TargetDivId;
        this.id = this.selection.m_tt_DivId;
        this.title = 'Range: ' + this.selection.getRangeStr();
        var range = this.selection.range;

        this.callParent(arguments);
        var iHeight = 22;
        
        var menuHandler = function(item) {
            this.hide();
            this.gview.pingClick('5-' + item.applog);
            switch (item.applog) {
                case 'Clip': this.seq2Clipboard(); break;
                case '1': this.modifyRange(); break;
                case '0': // zoom on range
                    this.gview.removeRangeSelection(range);
                    this.gview.startImageLoading(range[0], range[1] - range[0] + 1, {from_ui: true});
                    break;
                case '8': // zoom2seq
                    this.gview.removeRangeSelection(range);
                    this.gview.zoomSeq(Math.floor((range[1] + range[0])/2));
                    break;
                case '3': // marker
                    this.gview.m_App.addMarker(range);
                    this.gview.removeRangeSelection(range);
                    break;
                case '2': // new panel
                    this.gview.addNewPanel('', this.gview.getTotalSelectedRange());
                    this.highlightTargetRectangle();
                    break;
                case '4': this.gview.blastSelection(); break;
                case '5': this.gview.primerBlast(); break;
                case '6': // FASTA
                    this.gview.downloadData(true, "fasta", range); break;
                case '7': // genbank
                    this.gview.downloadData(true, "flat", range); break;
            }
        }
        var toolbar_items = [
            { text: 'Copy Sequence (Selection)', iconCls: 'xsv-search-loading', itemId: 'clip', disabled: true,
              tooltip: 'Copy selected sequence to clipboard', applog: 'Clip' },
            { text: 'Zoom On Range', iconCls: 'xsv-zoom_range', applog: '0' },
            { text: 'Zoom To Sequence', iconCls: 'xsv-zoom_seq', applog: '8' },
            { text: 'Set New Marker On Range', iconCls: 'xsv-markers', applog: '3' },
            { text: 'Modify Range', iconCls: 'xsv-zoom_range', applog: '1' }]

        if (this.gview.m_App.mf_MultiPanel) {
            toolbar_items.push({ text: 'Add New Panel On Range', iconCls: 'xsv-new_view', applog: '2'});
        }
        toolbar_items.push(
            { text: 'BLAST Search (Selection)', iconCls: 'xsv-blast', applog: '4' },
            { text: 'Primer BLAST (Selection)', iconCls: 'xsv-primer', applog: '5',
            disabled: (this.gview.m_App.m_ViewParams['acc_type'] !== 'DNA') },
            { text: 'Download FASTA (Selection)', iconCls: 'xsv-download', applog: '6' },
            { text: 'Download GenBank Flat File (Selection)', iconCls: 'xsv-download', applog: '7' }
        );
        if (this.gview.m_App.m_ViewParams.acc_type == 'protein') toolbar_items.pop();
        toolbar_items.forEach(function(itm, i) {
            Ext.apply(itm, {style: {margin: '0px',padding:'0px'}, x: 0, y: i * iHeight,
                scope: this, handler: menuHandler}); 
        }, this);

        var toolbar = new Ext.Toolbar({
            layout: 'absolute',
            border: false,
            items: toolbar_items,
            height: toolbar_items.length * iHeight, // + 16,
            width: 218
        });

        this.add(toolbar);
        if (this.gview.m_seqPrimers && range[1] - range[0] < 1024) this.getSelectedSeq([range[0], range[1]]);
        else this.down('#clip').setIconCls('xsv-clipboard');      
    },

    getSelectedSeq: function(range) {
        if (this.gview.getFlip()) range.push(range.shift());
        var self = this;
            app = this.gview.m_App;
        var callback = function(key) {
            var rq = new XMLHttpRequest();
            rq.open('get', app.m_CGIs.NetCache + key);
            rq.onload = function(e) {
                var item = self.down('#clip');
                item.enable();
                item.setIconCls('xsv-clipboard');
                var split = rq.response.split('\n');
                split.shift();
                self.tmpSelectedSeq = split.join('');
            }
            rq.send();
        }
        app.requestFile(this.gview, '?id=' + app.GI + '&ranges=' + range.join('-') + '&format=fasta', '', callback, true);
    },

    seq2Clipboard: function(){
        const el = document.createElement('textarea');
        el.textContent = this.tmpSelectedSeq;
        el.setAttribute('readonly', '');
        el.style.position = 'absolute';
        el.style.left = '-9999px';
        document.body.appendChild(el);
        el.select();
        document.execCommand('copy');
        document.body.removeChild(el);
    },

/*    seq2Clipboard: function(range) {
        if (this.gview.getFlip()) range.push(range.shift());
        var app = this.gview.m_App;
        var callback = function(key) {
            var rq = new XMLHttpRequest();
            rq.open('get', app.m_CGIs.NetCache + key);
            rq.onload = function(e) {
                var split = rq.response.split('\n');
                split.shift();
                navigator.permissions.query({name: "clipboard-write"}).then(function(result) {
                    console.log(result);
                    if (result.state == "granted" || result.state == "prompt") {
                        navigator.clipboard.writeText(split.join('')).then(function() {
                            console.log('Data copied');//
                        }, function() {
                            console.log('clipboard write failed');
                        });
                    } else {
                        const el = document.createElement('textarea');
                        el.textContent = split.join('');
                        el.setAttribute('readonly', '');
                        el.style.position = 'absolute';
                        el.style.left = '-9999px';
                        document.body.appendChild(el);
                        el.select();
                        document.execCommand('copy');
                        document.body.removeChild(el);
                    }
                });
            }
            rq.send();
        }
        app.requestFile(this.gview, '?id=' + app.GI + '&ranges=' + range.join('-') + '&format=fasta', '', callback);
    },*/

    modifyRange: function() {
        Ext.MessageBox.prompt('Modify Range', 'Enter New Range:', function(btn, text)
        {
            if (btn!='ok' || text.length==0) return;
            this.gview.m_App.handlePos(text, {
                allow_equal: true,
                ask_user: true,
                success: function(pos_range, options) {
                    if (pos_range.length !== 2) {
                        options.failure.call(options.scope, "Range required", options);
                        return;
                    }
                    this.selection.range[0] = pos_range[0];
                    this.selection.range[1] = pos_range[1];
                    this.selection.updateCoords();
                    this.setTitle('Range: ' + this.selection.getRangeStr());
                },
                failure: function(err_msg, options) {
                    Ext.MessageBox.alert('Error', err_msg);
                },
                scope: this
            });
        }, this, false, this.selection.getRangeStr(true));
    },

    onTargetOut: function(e){
        if (this.disabled || e.within(this.target.dom, true)){
            return;
        }
        this.highlightTargetRectangle();
        this.el.applyStyles('outline-style: none');
    },

    onTargetOver: function(e){
        this.highlightTargetRectangleMouseOver();
        if (this.pinned) {
            this.el.applyStyles('outline-style: solid;outline-color: red;outline-width: thin;');
        }
//        else this.show();
        this.show();
//        this.gview.ping({"jsevent":"mouseover","sv-event":"Range_ToolTip"});
    },

    onShow: function() {
        if (this.gview.m_ContextMenu && this.gview.m_ContextMenu.isVisible()) {
            this.gview.m_ContextMenu.hide();
        }
        this.callParent(arguments);
    },

    highlightTargetRectangleMouseOver: function() { this.target.dom.style.backgroundColor = '#94bfa2'; },
    highlightTargetRectangle: function() { this.target.dom.style.backgroundColor = '#4aa564'; }
});

//////////////////////////////////////////////////////////////////////////
SeqView.RangeSelection = function(view, range, event) {
    this.m_View = view;
    this.m_Resizing = false;
    if (!event) this.range = range;

    var parent_elem = Ext.fly(view.m_DivId);
    this.m_tt_TargetDivId = 'svtt-ruler-' + Ext.id();
    this.m_tt_DivId = 'svtt-' + Ext.id();
    var TargetDiv = '<div id=' + this.m_tt_TargetDivId + ' style="height:' + view.m_topRulerHeight + 'px; border: none; background-color: #95B9C7;"></div>';
    if (view.m_bottomRulerHeight)
        TargetDiv +='<div id=' + ('svtt-ruler-' + Ext.id()) + ' style="height:' + view.m_bottomRulerHeight + 'px; top:'
            + (view.m_Height - view.m_bottomRulerHeight - view.m_topRulerHeight)
            + 'px; position:relative; border:none; background-color: #95B9C7;"></div>';
    var RangeDiv = '<div class="range_selection sv-drag sv-highlight">'+ TargetDiv +'</div>';

    var tpl = new Ext.Template(RangeDiv);
    this.element = tpl.append(parent_elem, {}, true);
    this.m_ToolTip = null;

    this.element.on({ scope: this,
        'pointerdown': this.onPointerDown,
        'pointermove': this.onPointerMove,
        'click':     this.onClick
    });
    if (view.m_bottomRulerHeight) this.element.on({'mouseover': this.onMouseOver, scope: this});

    if (event /*means resizing*/) {
        var x = this.pageToViewX(Math.round(range[0]));
        this.range = [x, x];
        this.offsetX = 0;
        this.setElementPixRange(this.range);
        this.startResizing(event);
    } else {
        this.updateCoords();
    }
    if (!event && !this.m_ToolTip) {
        this.m_ToolTip = new SeqView.SelectedRangeToolTip({selection: this});
    }
};

SeqView.RangeSelection.prototype = {
    movePix: function(delta) { this.element.setLeft(this.element.getLeft(true) - delta); },
    pageToViewX: function(x) { return x - this.m_View.bodyDiv.getX(); },
    updateCoords: function() { this.setElementPixRange(this.rangeSeqToPix(this.range)); },

    getRangeStr: function(noLen) {
        var app = this.m_View.m_App;
        var start = app.posToLocalDisplay(this.range[0]);
        var end = app.posToLocalDisplay(this.range[1]);
        var len = noLen ? 0 : this.range[1] - this.range[0];
        return start + '..' + end + 
        (!len ? '' : 
             '<br><small>' + (len + 1).commify() + ' ' + app.m_ViewParams.units +'</small>');
    },

    remove: function() {
        this.element.remove();
        if (this.m_ToolTip) this.m_ToolTip.hide();
    },

    rangePixToSeq: function(range) {
        var view = this.m_View;
        var scp = view.m_ScrollPix;
        var over = view.getFlip() ? -0.5 : 0.5;
        var left = Math.round(view.pix2Seq(range[0]-scp, over));
        var right = Math.round(view.pix2Seq(range[1]-scp, -over));
        if (left <= right) 
            return [left, right];
        else
            return [right, left];
    },

    rangeSeqToPix: function(range) {
        var view = this.m_View;
        var left, right, over;
        if (!view.getFlip()) {
            over = 0.5;
            left = range[0];
            right = range[1];
        } else {
            over = - 0.5;
            left = range[1];
            right = range[0];
        }
        left = view.seq2PixScrolled(left - over);
        right = view.seq2PixScrolled(right + over);
        return [left, right];
    },

    setElementPixRange: function(range) {
        var view = this.m_View;
        // For logic underneath see markers.js:updateMarkerPos
        // TODO: move this code to view
        var left = Math.max(-10000, Math.min(10000, range[0]));
        var right = Math.max(-10000, Math.min(10000, range[1]));
        var width = right - left;
        if (width > 4000) {
            var mid_view = Math.round(view.getBodyWidth() / 2);
            if (Math.abs(left - mid_view) < Math.abs(right - mid_view)) {
                right = Math.max(mid_view + 2000, Math.min(right, left + 4000));
                left = Math.max(left, right - 4000)
            } else {
                left = Math.min(mid_view - 2000, Math.max(left, right - 4000));
                right = Math.min(right, left + 4000);
            }
            width = right - left;
        }
        this.element.setLeft(left).setWidth(width);
    },

    startResizing: function(e) {
        this.m_Resizing = true;
        this.m_View.m_App.eventHandlers('on', this);
    },

    onClick: function(e) {
        if (e.ctrlKey) {
            e.stopEvent();
            this.m_View.removeRangeSelection(true, this);
        }
    },

    onMouseOver: function(e) {
        var xy = e.getXY();
        if (this.m_ToolTip && e.target.id.indexOf('svtt-ruler') == 0) {
            if (this.m_ToolTip.target.id != e.target.id) {
                this.m_ToolTip.setTarget(e.target.id);
                this.m_ToolTip.show();
            }
        }
    },

    onPointerDown: function(e) {
        // if hit is close to the border, turn on resize
        this.time = e.time;
        var xy = e.getXY();
        var x = this.pageToViewX(xy[0]);
        var left = this.element.getLeft(true);
        var right = left + this.element.getWidth(true);
        var ld = Math.abs(x - left);
        var rd = Math.abs(x - right);
        if (ld < 5 || rd < 5) {
            e.stopPropagation();
            var pix_range = this.rangeSeqToPix(this.range);
            if (ld < rd) this.range = [pix_range[1], pix_range[0]];
            else this.range = pix_range;
            this.offsetX = x - this.range[1];
            this.m_View.popRangeElement(this);
            this.startResizing(e);
        }
    },

    onPointerUp: function(e) {
        if (e.time - this.time > 5000) {
            this.onClick(e);
            return;
        }
        var range = this.range;
        if (this.range[0] > this.range[1])
            range = [this.range[1], this.range[0]];
        
        this.range = this.rangePixToSeq([Math.max(0, range[0]), Math.min(range[1], this.m_View.m_Width)]);
        this.updateCoords();
        if (this.m_Resizing) {
            Ext.fly(this.m_View.m_DivId).setStyle('cursor', 'default');
            Ext.fly(e.getTarget()).setStyle('cursor', 'default');
        }

        this.m_Resizing = false;
        this.m_View.m_App.eventHandlers('un', this);

        // Let view decide whether it needs this range selection
        if (!this.m_View.addRangeSelection(this, this.range, range[1]-range[0])) return;

        if (!this.m_ToolTip)
            this.m_ToolTip = new SeqView.SelectedRangeToolTip({ selection: this, autoHide: true });
        if (this.m_ToolTip) {
            var el = Ext.get(this.m_tt_TargetDivId);
            if (el) {
                //this.m_ToolTip.show();
                this.m_ToolTip.showAt([e.clientX,e.clientY]);
                this.m_ToolTip.setTitle('Range: ' + this.getRangeStr());
                this.m_ToolTip.highlightTargetRectangleMouseOver();
            }
            this.m_View.m_App.ping({"jsevent":"click","sv-event":"Graph_Range_Selection"});
        }
    },

    onPointerMove: function(e) {
        NCBIGBUtils.ClearBrowserSelection();
        var x = this.pageToViewX(e.getX());
        x = Math.max(0, Math.min(x, this.m_View.getBodyWidth()));
        if (!this.m_Resizing) {
            var left = this.element.getLeft(true);
            var right = left + this.element.getWidth(true);
            var ld = Math.abs(x - left);
            var rd = Math.abs(x - right);
            if (ld < 5 || rd < 5) {
                Ext.fly(e.getTarget()).setStyle('cursor', 'e-resize');
            } else {
                Ext.fly(e.getTarget()).setStyle('cursor', 'default');
            }
            if (this.m_ToolTip && !this.m_ToolTip.isVisible()) {
                this.m_ToolTip.highlightTargetRectangle();
            }
            return;
        }
        e.ctrlKey = true; // forcing to show posTT
        this.m_View.updatePosTT(e);
        x -= this.offsetX;
        e.stopEvent();
        if (x == this.range[1]) return;
        this.range[1] = x;
        var left, width;
        if (this.range[0] < this.range[1]) {
            this.setElementPixRange(this.range);
        } else {
            this.setElementPixRange([this.range[1], this.range[0]]);
        }
    }
};
/*  $Id: markers.js 46963 2022-03-15 21:00:06Z borodine $
 * Authors:  Vlad Lebedev, Maxim Didenko, Victor Joukov
 * File Description: Visual markers
 *
 */

SeqView.MarkerFlags = {
    UserLock:      (1 << 0), ///< Locked by user
    SystemLock:    (1 << 1), ///< Locked by Sequence Viewer, unlockable by user
    Hollow:        (1 << 2), ///< Marker body is transparent
    PositionTitle: (1 << 3)  ///< Title is position-based and cleared if pos change
};
SeqView.MarkerNav = String.fromCharCode(171, 187);

//////////////////////////////////////////////////////////////////////
// SeqView.Marker

SeqView.Marker = (function () {
    var sm_Width = 8;

    // Coordinates are zero-based sequence position
    function constructor(minfo, range, name, flags, color) {
        this.m_MInfo = minfo;
        this.obj_coords = null;
        this.deleted = false;
        this.flags = flags;

        this.color = color || minfo.getNextColor();
        this.assignName(name ? name : minfo.suggestNextName());
        this.marker_num = minfo.m_NextMarkerNum++;

        this.span = 0;
        this.seq_pos = range[0];
        if (range[1])
            this.span = range[1] - range[0] + 1;

        this.marker_in_views = [];
        this.lines_in_views = [];
        this.addToAllViews();


    };

    constructor.getWidth = function() { return sm_Width; };

    var sm_PanoramaTmpl = new Ext.Template(
      '<div id="{id}" data-qtip="{qtip}" class="sv-marker sv-marker_base sv-opaque67">',
      '<div style="position:absolute;top:7px;left:0px;width:16px;height:2px;clip:rect(0 16px 2px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:4px;left:1px;width:14px;height:8px;clip:rect(0 14px 8px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:3px;left:2px;width:12px;height:10px;clip:rect(0 12px 10px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:2px;left:3px;width:10px;height:12px;clip:rect(0 10px 12px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:1px;left:4px;width:8px;height:14px;clip:rect(0 8px 14px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:0px;left:7px;width:2px;height:16px;clip:rect(0 2px 16px 0);background-color:{color};"></div>',
      '<div id="{id}_line" class="marker_line sv-opaque5" style="width:1px;border-left:1px solid {color};">',
          '<div id="{id}_body" class="sv-opaque3 sv-drag sv-highlight" style="background-color:{color};position:absolute;width:100%;height:100%"></div>',
      '</div>',
      '<div data-qtip="Locked" id="the_{id}_lock" style="display:{lock};" class="marker_locked_ov"></div></div>'
    );
    
    var sm_PanoramaRectTmpl = new Ext.Template(
      '<div id="{id}" data-qtip="{qtip}" class="sv-marker sv-marker_base sv-opaque67">',
      '<div style="position:absolute;top:-2px;left:2px;width:13px;height:14px;background-color:{color};opacity:0.4;"></div>',
      '<div id="{id}_line" class="marker_line sv-opaque5" style="top:14px;width:1px;border-left:1px solid {color};">',
      '</div></div>'
    );
    

    var sm_GraphicTmpl = new Ext.Template(
      '<div id="{id}" class="sv-marker sv-marker_half_base" style="top:0px;">',
      '<div style="position:absolute;top:7px;left:0px;width:8px;height:2px;clip:rect(0 8px 2px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:4px;left:1px;width:7px;height:8px;clip:rect(0 7px 8px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:3px;left:2px;width:6px;height:10px;clip:rect(0 6px 10px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:2px;left:3px;width:5px;height:12px;clip:rect(0 5px 12px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:1px;left:4px;width:4px;height:14px;clip:rect(0 4px 14px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:0px;left:7px;width:1px;height:16px;clip:rect(0 1px 16px 0);background-color:{color};"></div>',
      '<div id="{id}_line" class="marker_line sv-opaque5" style="width:0px;border-left:1px solid {color};border-right:1px solid {color};">',
          '<div id="{id}_body" class="sv-opaque3 sv-drag sv-highlight" style="background-color:{color};position:absolute;width:100%;height:100%"></div>',
      '</div>',
      '<div id="the_{id}_back" class="marker_label_back sv-opaque67" style="background-color:{color};width:{label_length}px;"></div>',
       '<div id="{id}_label" class="marker_label">{trimmed_label}',
      ' <img aria-label="icon" id="the_{id}_lock" style="margin-top:1px;display:{lock};" src="' + SeqView.base_url + 'images/lock.png"></div></div>'
    );
    
    var sm_GraphicRectTmpl = new Ext.Template(
      '<div id="{id}" class="sv-marker sv-marker_half_base" style="top:0px;">',
      '<div id="{id}_line" class="marker_line sv-opaque5" style="width:0px;border-left:1px solid {color};border-right:1px solid {color};">',
          '<div id="{id}_body" class="sv-opaque3 sv-drag sv-highlight" style="background-color:{color};position:absolute;width:100%;height:100%"></div>',
      '</div>',
      '<div id="the_{id}_back" class="marker_label_back sv-opaque67" style="background-color:{color};opacity:0.4;top:1px;left:1px;width:18px;"></div>',
      '<div id="{id}_label" class="marker_label" style="font-size:13px;left:1px;top:0px;">' + SeqView.MarkerNav + '</div>',
      '</div>'
    );


    var sm_GraphicHollowTmpl = new Ext.Template(
      '<div id="{id}" class="sv-marker sv-marker_half_base" style="top:0px;">',
      '<div style="position:absolute;top:7px;left:0px;width:8px;height:2px;clip:rect(0 8px 2px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:4px;left:1px;width:7px;height:8px;clip:rect(0 7px 8px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:3px;left:2px;width:6px;height:10px;clip:rect(0 6px 10px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:2px;left:3px;width:5px;height:12px;clip:rect(0 5px 12px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:1px;left:4px;width:4px;height:14px;clip:rect(0 4px 14px 0);background-color:{color};"></div>',
      '<div style="position:absolute;top:0px;left:7px;width:1px;height:16px;clip:rect(0 1px 16px 0);background-color:{color};"></div>',
      '<div id="{id}_line" class="marker_line sv-opaque5" style="width:0px;border-top:5px solid {color};border-left:1px solid {color};border-right:1px solid {color};">',
          '<div id="{id}_body" class="sv-transparent sv-drag sv-highlight" style="background-color:{color};position:absolute;width:100%;height:100%"></div>',
      '</div>',
      '<div id="the_{id}_back" class="marker_label_back sv-opaque67" style="background-color:{color};width:{label_length}px;"></div>',
      '<div id="{id}_label" class="marker_label">{trimmed_label}',
      ' <img aria-label="icon" id="the_{id}_lock" style="margin-top:1px;display:{lock};" src="' + SeqView.base_url + 'images/lock.png"></div></div>'
    );

    constructor.getTemplate = function(is_panorama, fill_bg, rect) {
        if (is_panorama) return rect ? sm_PanoramaRectTmpl : sm_PanoramaTmpl;
        if (rect) return sm_GraphicRectTmpl;
        return fill_bg ? sm_GraphicTmpl : sm_GraphicHollowTmpl;
    };

    return constructor;
}) ();


SeqView.Marker.prototype = {

    displayPos: function(forPDF) {
        var pos = forPDF ? this.seq_pos : this.seq_pos + 1;
        if (this.span) pos += ':' + (this.seq_pos + this.span);
        return pos;
    },

    displayLocalPos: function() {
        var app = this.m_MInfo.m_App;
        var pos = app.posToLocal(this.seq_pos);
        if (this.span) {
            if (app.getFlip()) {
                pos = app.posToLocalDisplay(this.seq_pos + this.span - 1) + ':' + pos;
            } else {
                pos += ':' + app.posToLocalDisplay(this.seq_pos + this.span - 1);
            }
        }
        return pos;
    },

    addToAllViews: function() {
        var app = this.m_MInfo.m_App;
        app.forEachView(function(view) { this.addToView(view); }, this);
        if (app.m_TextView)
            app.m_TextView.AddOrUpdateMarker(this); // text view
        app.fireEvent('marker_created', this);
    },

    getCreateParams: function(is_panorama, idx) {
        var elem_id = 'marker_' + idx + '_' + this.marker_num;
        var user_lock = this.flags & SeqView.MarkerFlags.UserLock;
        this.nvMarker = this.marker_name == SeqView.MarkerNav;
        return {
            template: SeqView.Marker.getTemplate(is_panorama,
                !(this.flags & SeqView.MarkerFlags.Hollow), this.nvMarker),
            options: {
                id: elem_id, color: this.color,
                num: this.marker_num, trimmed_label: this.marker_name_trimmed,
                label_length: user_lock ? this.marker_name_back_lock : this.marker_name_back,
                lock: user_lock ? 'inline' : 'none', prefix: this.m_MInfo.m_App.m_CGIs.baseURL
            }
        };
    },

    addToView: function(view) {
        if (view.isAlignment())
            return;

        var elem_id = 'marker_' + view.m_Idx + '_' + this.marker_num;

        var marker = view.createMarkerElem(this);
        this.m_MInfo.m_AllMarkers[elem_id] = this;

        var marker_line = Ext.get(elem_id + '_line');
        marker_line.setHeight(view.getHeight(true));

        this.lines_in_views[view.m_Idx] = marker_line;
        this.marker_in_views[view.m_Idx] = marker;

        this.updateMarkerPos(view);

        var marker_label = Ext.get(elem_id + '_label');
        if (marker_label) {
            marker_label.on({
                'pointerdown': this.onPointerDown,
                'pointermove': this.onPointerMove,
                'contextmenu': this.onContextMenu,
                scope: this
            });
        }
        //creating tooltip object that targets this marker
        if (view.m_Type == "graphical") {
            this.m_ToolTip = new SeqView.MarkerToolTip({
                selection: this,
                target: marker_label,
                gview: view
            });
        }

    },
    removeFromView: function(view) {
        //removing marker tooltips
        if (this.m_ToolTip) {
//            this.m_ToolTip.remove();
            this.m_ToolTip.destroy();
            delete this.m_ToolTip;
        }
        var elem_id = 'marker_' + view.m_Idx + '_' + this.marker_num;
        if (this.marker_in_views[view.m_Idx]) {
            this.marker_in_views[view.m_Idx].remove();
            delete this.m_MInfo.m_AllMarkers[elem_id];
        }
    },

    lockMarker: function() {
        if (this.flags & SeqView.MarkerFlags.SystemLock) return;
        this.flags = this.flags ^ SeqView.MarkerFlags.UserLock;
        if (!(this.flags & SeqView.MarkerFlags.UserLock) &&
             (this.flags & SeqView.MarkerFlags.PositionTitle)) {
            this.flags = this.flags ^ SeqView.MarkerFlags.PositionTitle;
            this.rename(this.m_MInfo.suggestNextName());
        }

        // update views
        this.m_MInfo.m_App.forEachView(function(view) {
            var elem_id = 'the_marker_' + view.m_Idx + '_' + this.marker_num + '_lock';
            var back_id = 'the_marker_' + view.m_Idx + '_' + this.marker_num + '_back';
            var the_marker = Ext.get(elem_id);
            var the_back = Ext.get(back_id);

            // it's a bit tricky. Marker backgroung is transparent and the label is opaque
            // this means we have to use absolute positioning and resize the background manually to put/remove the lock icon
            if (the_marker) {
                var cur_display = the_marker.getStyle('display');
                the_marker.setStyle('display', cur_display == 'none' ? 'inline' : 'none');
                if (!view.isPanorama()) {
                    the_back.setStyle('width', '' + ((this.flags & SeqView.MarkerFlags.UserLock) ? this.marker_name_back_lock : this.marker_name_back) + 'px');
                }
            }
        }, this);

        // Text View
        if (this.m_MInfo.m_App.m_TextView)
            this.m_MInfo.m_App.m_TextView.UpdateLockMarker(this);
    },

    deleteMarker: function() {
        this.m_MInfo.m_App.forEachView(function(view) {
            this.removeFromView(view);
        }, this);

        // remove from text view
        if (this.m_MInfo.m_App.m_TextView)
            this.m_MInfo.m_App.m_TextView.RemoveMarker(this);

        this.deleted = true;
        this.m_MInfo.m_App.fireEvent('marker_deleted', this);

    },


    getMarker: function(elem_id) {
        var s = elem_id.split('_');
        view = this.m_MInfo.m_App.findView(s[1]);
        marker = this.marker_in_views[s[1]];
        return { view: view, marker: marker };
    },

    setCursor: function(elem_id, cursor) {
        this.getMarker(elem_id).marker.setStyle('cursor', cursor);
    },

    getSeqPos: function(elem_id, pix_pos) {
        var s = elem_id.split('_');
        var view = this.m_MInfo.m_App.findView(s[1]);
        if (view) {
            var pix = pix_pos - view.m_ScrollPix;
            return view.pix2Seq(pix);
        }
        return 0;
    },

    // set marker name and calculate derived names and parameters
    assignName: function(name) {
        this.marker_name_trimmed = name.trimToPix(100);
        this.marker_name = name;
        this.visible_name = NCBIGBUtils.sanitize(name);
        this.marker_name_back = this.marker_name_trimmed.visualLength() + 2;
        this.marker_name_back_lock = this.marker_name_back + 16;
    },

    setName: function(name) {
        this.assignName(name);

        this.m_MInfo.m_App.forEachView(function(view) {
            if (view.isAlignment())
                return;

            var elem_id = 'marker_' + view.m_Idx + '_' + this.marker_num;
            var back_id = 'the_marker_' + view.m_Idx + '_' + this.marker_num + '_back';
            var elem_label_id = 'marker_' + view.m_Idx + '_' + this.marker_num + '_label';

            var marker_div = Ext.get(elem_id);
            var the_back = Ext.get(back_id);
            var marker_label_div = Ext.get(elem_label_id);
            if (marker_label_div) {// no labels in overview
                var lock_html = marker_label_div.dom.innerHTML;
                lock_html = lock_html.slice(lock_html.indexOf('<'), lock_html.length);

                marker_label_div.update(this.marker_name_trimmed + ' ' + lock_html);
                the_back.setStyle('width', '' + ((this.flags & SeqView.MarkerFlags.UserLock) ? this.marker_name_back_lock : this.marker_name_back) + 'px');
            }
        }, this);
        var marker_seq = Ext.get('seqtext-marker_' + this.m_MInfo.m_App.m_Idx + '_' + this.marker_num);
    },

    rename: function(name) {
        this.setName(name);
        //changing tooltip title too
        var title = '&nbsp;&nbsp;&nbsp;' + name;
        this.m_ToolTip.setTitle(title);
        //
        this.m_MInfo.updateInfo();
    },

    setSeqPos: function(pos, reload) {
        this.seq_pos = pos;

        this.m_MInfo.m_App.forEachView(function(view) {
            if (view.isAlignment())
                return;
            this.updateMarkerPos(view);
        }, this);

        if (reload) {
            this.m_MInfo.updateInfo();
        }

        if (reload && this.m_MInfo.m_App.m_TextView) { // text view
            this.m_MInfo.m_App.m_TextView.SetSeqPos(this);
        }
    },

    setSeqRange: function(range, reload) {
        var pos = range[0];
        if (range[1] && range[1] > pos)
            this.span = range[1] - pos + 1; // inclusive range
        else
            this.span = 0;
        this.setSeqPos(pos, reload);
    },

    setTop: function(h) { this.marker.setTop(h); },

    movePix: function(elem_id, delta, validate) {
        var vm = this.getMarker(elem_id);
        var new_left = vm.marker.getLeft(true) - delta;
        var marker_width = SeqView.Marker.getWidth();
        if (validate) { // validate only when dragging with mouse
            new_left = Math.max(-marker_width, new_left);
            var rightmost_pix = vm.view.getMostRightPix();
            new_left = Math.min(rightmost_pix /* - vm.marker.getWidth() */ + marker_width, new_left);
        }
        var pix = new_left + marker_width - vm.view.m_ScrollPix;
        var flip = vm.view.getFlip && vm.view.getFlip(); // ugly, but overview does not have a getFlip
        var new_seq_pos = vm.view.pix2Seq(pix) - (flip ? this.span : 0);
        // Check new position in sequence space and fix it if possible
        var fixed_new_seq_pos = Math.max(0, new_seq_pos);
        fixed_new_seq_pos = Math.min(vm.view.m_App.m_SeqLength - this.span - 1, fixed_new_seq_pos);
        this.seq_pos = fixed_new_seq_pos;
        if (fixed_new_seq_pos !== new_seq_pos) {
            new_left = vm.view.seq2Pix(fixed_new_seq_pos);
        }
        vm.marker.setLeft(new_left);

        if (validate) { // syncronize other markers
            this.m_MInfo.m_App.forEachView(function(view) {
                if (view.isAlignment() || view.m_Idx == vm.view.m_Idx)
                    return;
                this.updateMarkerPos(view);
            }, this);

            if (this.m_MInfo.m_App.m_TextView) { // text view
                this.m_MInfo.m_App.m_TextView.AddOrUpdateMarker(this);
            }
        }
    },

    scrollPix: function(view, delta) {
        if (view.isAlignment())
            return;

        var marker = this.marker_in_views[view.m_Idx];
        var new_left = marker.getLeft(true) - delta;
        marker.setLeft(new_left);
    },

    updateMarkerSize: function(view) {
        if (view.isAlignment())
            return;
        var idx = view.m_Idx;
        var marker_line;
        if (this.lines_in_views[idx] != null) {
            marker_line = this.lines_in_views[idx];
        } else {
            this.addToView(view);
            marker_line = this.lines_in_views[idx];
        }
        marker_line.setHeight((view.m_FromCgi.img_height || view.m_Height) - 18);
    },

    updateMarkerPos: function(view) {
        if (view.isAlignment()) return;
        var view_idx = view.m_Idx;
        var marker = this.marker_in_views[view_idx];
        if (!marker) return;
        var effective_pos = this.seq_pos - 0.5;
        var pix_pos = view.seq2PixScrolled(effective_pos);
        var marker_line = this.lines_in_views[view_idx];
        var end_pos = view.seq2PixScrolled((this.span || 1)  + effective_pos);
        if (end_pos < pix_pos) {
            var pos = end_pos;
            end_pos = pix_pos;
            pix_pos = pos;
        }
        var pix_width = end_pos - pix_pos;
        // small correction for visible offset
        pix_pos -= 1;
        end_pos += 1;
        if (!this.span && pix_width < 9) {
            pix_pos = view.seq2PixScrolled(this.seq_pos);
            pix_width = 1;
            marker_line.setStyle("border-right", "none");
        }
        else marker_line.setStyle("border-right", "1px solid " + this.color);
        // Some browsers, WebKit-based in particular, don't like coordinates too far
        // beyond the window. We need to have an ability to hide a marker, that is why we
        // we cull coordinate by +- 10000
        pix_pos = Math.max(-10000, Math.min(10000, pix_pos - SeqView.Marker.getWidth()));
        if (pix_width > 1) {
            end_pos = Math.max(-10000, Math.min(10000, end_pos - SeqView.Marker.getWidth()));
            // Also, we should observe IE limit for width of transparent objects, ca. 4100px
            // Being on the safe side, also compatible with out image size - 4000px
            // Pick the end closest to the mid-window and adjust the other end so that
            // total width is less than 4000.
            if (end_pos - pix_pos > 4000) {
                var mid_view = Math.round(view.getBodyWidth() / 2);
                if (Math.abs(pix_pos - mid_view) < Math.abs(end_pos - mid_view)) {
                    end_pos = Math.max(mid_view + 2000, Math.min(end_pos, pix_pos + 4000));
                    pix_pos = Math.max(pix_pos, end_pos - 4000)
                } else {
                    pix_pos = Math.min(mid_view - 2000, Math.max(pix_pos, end_pos - 4000));
                    end_pos = Math.min(end_pos, pix_pos + 4000);
                }
            }
            pix_width = end_pos - pix_pos;
        }
        marker.setLeft(pix_pos);
        marker_line.setWidth(pix_width);
        // To satisfy IE pre-8
        if (pix_width > 1) {
            marker.setWidth(pix_width + 16);
        }
    },


    onPointerDown: function(e) {
        if (e.pointerType == 'mouse' && e.button !== 0) return;
        if ((this.flags & (SeqView.MarkerFlags.SystemLock | SeqView.MarkerFlags.UserLock)) || Ext.fly(e.getTarget()).hasCls('sv-drag')) {
            delete this.offsetX;
            // ignore event for locked marker and pass it to the other handlers
            return;
        }
        //hide tooltip when drugging marker
        if (this.m_ToolTip && !this.m_ToolTip.isPinned()) this.m_ToolTip.hide();
        var marker_node = Ext.fly(e.getTarget()).findParent('div.sv-marker');
        if (!marker_node) return;

        this.curVM = this.getMarker(e.getTarget().id);
        this.offsetX = this.curVM.marker.getLeft(true) - e.getX();
        this.curVM.marker.setStyle('cursor', 'move');
        this.m_MInfo.m_App.eventHandlers('on', this);
        if (Ext.isIE) e.stopPropagation(); else e.stopEvent();
    },

    onPointerUp: function(e) {
        e.stopEvent();
        if (this.curVM) this.curVM.marker.setStyle('cursor', 'pointer');
        delete this.curVM;
        this.m_MInfo.m_App.eventHandlers('un', this);
        this.m_MInfo.updateInfo();
        //updating marker tooltip content
        this.m_ToolTip.updateMarkerToolTipContent(this.span, this.seq_pos);
    },

    onPointerMove: function(e) {
        if (!this.curVM) return; 
        NCBIGBUtils.ClearBrowserSelection();
        var vm = this.curVM;

        var marker_width = SeqView.Marker.getWidth();
        var flip = vm.view.getFlip && vm.view.getFlip(); // ugly, but overview does not have a getFlip
        // validate coordinate against visible window
        var new_left = Math.max(-marker_width, e.getX() + this.offsetX);
        var rightmost_pix = vm.view.getMostRightPix();
        var dWidth = vm.marker.getWidth() - marker_width;
        new_left = Math.min(rightmost_pix - dWidth, new_left);
        var pix = new_left + marker_width - vm.view.m_ScrollPix + ((dWidth && !this.span) ? vm.view.m_bpPix/2 : 0);
        var new_seq_pos = vm.view.pix2Seq(pix) - (flip ? this.span : 0);
        // Check new position in sequence space and fix it if possible:
        // taking into account vm.marker.getWidth() (above) takes care of this
        // approximately. Should work worse for longer sequences where
        // ration between sequence space and pixel space is higher.
        // So we adjust new_seq_pos so that both ends of marker are in
        // valid range.
        var fixed_seq_pos = Math.max(0, new_seq_pos);
        fixed_seq_pos = Math.min(vm.view.m_App.m_SeqLength - this.span, fixed_seq_pos);
        this.seq_pos = fixed_seq_pos;
        vm.marker.setLeft(new_left);
        // syncronize other markers
        this.m_MInfo.m_App.forEachView(function(view) {
            if (view.isAlignment() || view.m_Idx == vm.view.m_Idx)
                return;
            this.updateMarkerPos(view);
        }, this);

        if (this.m_MInfo.m_App.m_TextView) { // text view
            this.m_MInfo.m_App.m_TextView.AddOrUpdateMarker(this);
        }
    },


    onContextMenu: function(e) {
        e.stopEvent();
    },

    //////////////////////////////////////////////////////////////////////////
    // showSequence:

    showSequence: function() {
        this.m_MInfo.m_App.createTextView(this.seq_pos);
    },


    centerInView: function(view) {
        var new_len = view.m_VisLenSeq;
        var new_from = this.seq_pos - view.m_VisLenSeq / 2;
        // Center range of marker, not its beginning if possible
        if (this.span) {
            if (this.span < new_len)
                new_from += Math.floor(this.span / 2);
            else
                new_from += view.m_VisLenSeq / 3;
        }

        if (new_from < 0) new_from = 0;
        if (new_from + new_len > view.m_App.m_SeqLength) new_len = view.m_App.m_SeqLength - new_from;
        if (new_from === view.m_VisFromSeq && new_len === view.m_VisLenSeq) return;
        //view.startImageLoading(new_from, new_len);
        view.moveTo(new_from, new_len);
    },


    zoomSeqMarker: function(view) { view.zoomSeq(this.seq_pos); }
};


/********************************************************************/
//////////////////////////////////////////////////////////////////////
// SeqView.MarkersInfo
/********************************************************************/


SeqView.MarkersInfo = (function () {

    function constructor(app) {
        this.m_App = app;
        this.m_Dlg = null;
        this.m_AllMarkers = {};
        this.m_AllMarkersRefs = [];
        this.m_MarkerColorIdx = 0;
        this.m_NextMarkerNum = 1;
    };

    constructor.deleteMarker = function(app_idx, marker_num) {
        var app = SeqView.App.findAppByIndex(app_idx);
        if (app && app.m_MarkersInfo) {
             app.m_MarkersInfo.deleteMarker(marker_num);
        }
    };

    constructor.setMarkerPosition = function(app_idx, marker_num) {
        var app = SeqView.App.findAppByIndex(app_idx);
        if (app && app.m_MarkersInfo) {
             app.m_MarkersInfo.setMarkerPosition(marker_num);
        }
    };

    constructor.showMarker = function(app_idx, marker_num, zoom) {
        var app = SeqView.App.findAppByIndex(app_idx);
        if (app && app.m_MarkersInfo) {
             var marker = app.m_MarkersInfo.getMarker(marker_num);
             if (marker) app.forEachView(function(view) {
                 if (view.isGraphic()) marker[zoom ? 'zoomSeqMarker' : 'centerInView'](view);
             });
        }
    };

    constructor.renameMarker = function(app_idx, marker_num) {
        var app = SeqView.App.findAppByIndex(app_idx);
        if (app && app.m_MarkersInfo) {
             app.m_MarkersInfo.renameMarker(marker_num);
        }
    };

    return constructor;
}) ();

SeqView.MarkersInfo.prototype = {
    getNextColor: function() {
        var sm_MarkerColors = ['blue', 'red', 'green', 'cyan', '#FF00FF', 'black'];
        var color = sm_MarkerColors[this.m_MarkerColorIdx++];
        if (this.m_MarkerColorIdx > sm_MarkerColors.length-1)
            this.m_MarkerColorIdx = 0;
        return color;
    },

    suggestNextName: function() {
        //var name = this.m_NextMarkerNum;
        var array = this.m_AllMarkersRefs;
        var marker_nums = [];
        for (var i = 0, len = array.length; i < len; i++) {
            var m = array[i];
            if (!m.deleted) {
                var parts = m.marker_name.match(/^Marker (\d+)$/);
                if (parts != null) {
                    marker_nums.push(parts[1]);
                }
            }
        }
        marker_nums.sort(function(a,b){return a-b});
        var next_num = 1;
        for (var i = 0, len = marker_nums.length; i < len; i++) {
            if (next_num != marker_nums[i]) break;
            next_num++;
        }
        return "Marker " + next_num;
 },

//////////////////////////////////////////////////////////////////////////
// showDlg:

    showDlg: function(view) {
        var active_marker_count = 0;
        this.forEachMarker(function() { active_marker_count++; });

        if (active_marker_count == 0) {
            this.MarkerDlg(view);
            return;
        }
        this.m_Dlg = new Ext.Window({
            layout:'fit',
            title:'Markers',
            minWidth:480, width:600, height:350,
            cls: 'SeqViewerApp',
            closeAction: 'destroy',
            constrain: true, 
            collapsible: true,
            listeners: {scope: this,
                close: function() {
                    delete this.m_Dlg;
                    this.m_App.un('origin_changed', this.onOriginChanged, this);
                    this.m_App.un('strand_changed', this.onStrandChanged, this);
                    this.m_App.resizeIFrame();
                    this.m_App.m_DialogShown = false; 
                }
            },
            tbar: [
                {text: 'Add Marker', iconCls:'xsv-marker-add', tooltip:'Add New Marker', scope:this,
                 handler: function() { this.MarkerDlg(view); }},
                '->',
                {text:'Remove all markers', iconCls:'xsv-clear_markers', scope:this,
                 handler: function() {
                     Ext.MessageBox.confirm( 'Confirm', 'Clear all markers?', function(btn) {
                         if (btn == 'yes') this.reset();
                     }, this);
                 }}
            ],
            items: [{xtype:'panel', autoScroll: true, itemId: 'infoPanel', layout:'fit'}],
            buttons: [
                {text: 'Download', scope: this, handler: function(button) { view.saveFile(this.getMarkerInfoURL(true), '7-2-1'); }},
                {text: 'Close', scope: this, handler: function() { this.m_Dlg.close(); }}
            ]
        });
        this.m_App.on('origin_changed', this.onOriginChanged, this);
        this.m_App.on('strand_changed', this.onStrandChanged, this);
        this.m_App.resizeIFrame(400);
        this.m_App.m_DialogShown = true;
        this.updateInfo();    
        this.m_Dlg.show();
    },

    onOriginChanged: function(app) { this.updateInfo(); },
    onStrandChanged: function(app) { if (app.getOrigin() != 0) this.updateInfo(); },

    addMarker: function(range, name, lock, color, flags) {
        var newMarker = new SeqView.Marker(this, range, name, lock, color, flags);
        this.m_AllMarkersRefs.push(newMarker);
        this.updateInfo();
        return newMarker;
    },

//////////////////////////////////////////////////////////////////////////
// addMarkerByData:

    addMarkerByData: function(data) {
        this.addMarker(data[0], data[1], data[2], data[3], data[4]);
    },

    MarkerDlg: function(view, pos, oldMarker){
        var setMarkerDlg;

        if (view) { // pos is X coordinate on the image
            pos = pos || Math.floor(view.getBodyWidth()/2);
            pos = this.m_App.posToLocalDisplay(view.pix2Seq(pos - view.m_ScrollPix));
        } else { // if pos set it's position on the sequence
            pos = pos || this.m_App.posToLocalDisplay(Math.floor(this.m_App.m_SeqLength / 2));
            view = this;
        }

        setMarkerDlg = new Ext.Window({
            modal:true, title:'Set New Marker',
            constrain: true,
            app: view.m_App,
            width:270,
            collapsible:false, resizable:false, closeAction:'close', iconCls:'xsv-markers',
            itemId: 'setMarkerDlg',
            border: false,
            buttonAlign:'center',
            viewModel: { data: {
                color: oldMarker ? oldMarker.color : this.getNextColor(),
                lock: false,
                name: oldMarker ? oldMarker.visible_name : this.suggestNextName(),
                position: pos
            }},
            items:[
               {xtype: 'form', frame: true,
                items:[
                   {xtype:'textfield', width:'100%', fieldLabel: 'Name', hidden: oldMarker != null, bind: '{name}'},
                   {xtype:'textfield', width:'100%', fieldLabel: 'Pos/Range',  bind: '{position}'},
                   {xtype: 'container', layout: 'hbox',  height: 22,
                    items: [
                        {xtype:'displayfield', fieldLabel: 'Color'},
                        {xtype:'component', width: 22, height: 22, bind: {style: {backgroundColor:  '{color}'}}},
                        {xtype:'button',
                            menu: new Ext.menu.ColorPicker(
                              {listeners: { 'select': function(picker, choice) {
                                    picker.up('#setMarkerDlg').getViewModel().setData({color: '#' + choice});
                                }}})}]
                   },
                   {xtype:'checkbox', height:22, labelSeparator:'', boxLabel:'Lock Marker', hidden: oldMarker != null, bind: '{lock}'}
                ]}
            ],
            buttons:[
               {text:'OK', scope:this, handler: function() {
                  var data = setMarkerDlg.getViewModel().getData();
                  if (!data.name ||  data.name.length > 50) {
                      Ext.MessageBox.show({ title: 'Set New Marker', msg: 'Invalid marker name.',
                                            buttons: Ext.MessageBox.OK, icon:Ext.MessageBox.ERROR});
                      return;
                  }
                  this.m_App.handlePos(data.position.toString(), {
                      allow_equal: true,
                      ask_user: true,
                      success: function(pos_range, options) {
                          setMarkerDlg.close();
                          if (oldMarker) this.deleteMarker(oldMarker.marker_num);
                          this.addMarker(pos_range, data.name, data.lock, data.color);
                      },
                      failure: function(err_msg, options) {
                          Ext.MessageBox.show({title: 'Set New Marker', msg: err_msg,
                              buttons: Ext.MessageBox.OK, icon: Ext.MessageBox.ERROR });
                      },
                      scope: this
                  });
               }},
               {text:'Cancel', handler: function() {setMarkerDlg.close(); } }
             ]
        });
        setMarkerDlg.show();
    },


    getMarkerInfoURL: function(for_download, format) {
        var url = this.m_App.m_CGIs.ObjCoords + '?objcoords=1&id=' + this.m_App.GI 
        var positions = [];
        var names = [];
        this.forEachMarker(function(m) {
            positions.push(m.seq_pos);
            names.push(m.visible_name);
            if (m.span) {
                positions.push(m.seq_pos + m.span - 1);
                names.push(m.visible_name);
            }
        });
        if (for_download) {
            url += '&download=true&fmt=' + (format ? format : 'csv');
            for (var i = 0; i < positions.length; i++) {
                url += '&pos=' + positions[i] + '&name=' + names[i];
            }
        } else {
            var sorted = positions.sort(function(a,b) { return a - b;});
            for (var i = 0, prev = sorted[0] - 1; i < sorted.length; i++) {
                if (sorted[i] != prev) url += '&pos=' + sorted[i];
                prev = sorted[i];
            }
        }
        return url;
    },


    parseObjCoords: function(obj_coords_json, marker) {
        var data = obj_coords_json,
            pos  = marker.seq_pos,
            end_pos,
            fnParams = this.m_App.m_Idx + ',' + marker.marker_num;
        if (marker.span)
            end_pos = pos + marker.span - 1;
        var html = '';
        var useHGVS = this.m_App.m_ViewParams.organism == "Homo sapiens";
        var cur_pos = -1;
        for (var i = 0; i != data.length; i++) {
            var row_data = data[i];
            var marker_pos = row_data['marker_pos'];
            if (marker_pos != pos  &&  marker_pos != end_pos) continue;
            if (html.length == 0) {
                html += '<table class="xsv-markerinfo" width="98%">'
                html += '<tr><th><b>Name:</b> ' + marker.visible_name;
                if (!marker.nvMarker) html += '<br><a href="#" onClick="SeqView.MarkersInfo.renameMarker(' + fnParams + ');">Edit</a></th>';
                html += '<th width="75%" colspan=' + (useHGVS ? 4 : 3) + '>' + marker.m_ToolTip.rangeText.split(' ', 1)[0]
                    + '<a href="#" onClick="SeqView.MarkersInfo.showMarker(' + fnParams + ', true);" title="Zoom To Sequence At Marker">'
                    + marker.m_ToolTip.rangeText.split(':')[1] + '</a>';
                if (!marker.nvMarker) html += ' <br>' + //<a href="#" onClick="SeqView.MarkersInfo.showMarker(' + fnParams + ', true);" style="float:left;">&nbsp;ZoomTo</a>' +
                    '<a href="#" onClick="SeqView.MarkersInfo.setMarkerPosition(' + fnParams + ');">Edit</a>';
                html += '<a href="#" onClick="SeqView.MarkersInfo.deleteMarker(' + fnParams + ');" class="xsv-marker_remove_link">Remove</a>';
                html += '</th></tr>';
                if (useHGVS)
                    html +='<tr><th width="25%">Accession/Locus tag</th><th width="10%">Location</th><th width="10%">Relative to</th><th width="30%"><a href="http://varnomen.hgvs.org" target="_blank">HGVS Name</a></th><th width="25%">Sequence</th></tr>';
                else
                    html +='<tr><th width="25%">Accession/Locus tag</th><th width="20%">Location</th><th width="10%">Relative to</th><th width="45%">Sequence</th></tr>';
            }
            html += this.addMarkersRow(row_data, useHGVS, false, cur_pos != marker_pos);
            if (cur_pos != marker_pos  &&  this.m_App.m_Origin)
                html += this.addMarkersRow(row_data, useHGVS, true);
            cur_pos = marker_pos;
        }
        html += '</table>';
        return html;
    },

    addMarkersRow: function(data, useHGVS, use_local, use_separator) {
        var pos_mapped = '';
        if (typeof(data['pos_mapped']) != 'undefined') {
            pos_mapped = data['pos_mapped'];
            if (use_local) {
                pos_mapped = this.m_App.posToLocal(pos_mapped);
            } else {
                if (pos_mapped >= 0) pos_mapped += 1;
            }
        }
        if (!useHGVS && pos_mapped.length == 0) return '';

        var html = use_separator ? '<tr class="sv-rowgroup">' : '<tr>';
        html += '<td>' + data['title'] + '</td>'; // Accession
        html += '<td>' + pos_mapped + '</td>'; // Location
        var hgvs = data['hgvs_position'];
        var rel_to = "Seq start";
        if (use_local) {
            rel_to = "Current origin";
        }
        else if (hgvs  &&  hgvs.charAt(0) == 'c') {
            rel_to = "CDS start";
        }
        html += '<td>' + rel_to + '</td>'; // Relative to
        if (useHGVS) {
            if (hgvs  &&  data['title']) {
                hgvs = data['title'] + ':' + hgvs;
            }
            html += '<td>' + hgvs + '</td>';
        }
        html += '<td><span style="font-family:monospace">' + data['sequence']+'</span></td>';
        html += '</tr>';
        return html;
    },


    updateInfo: function() {
        if (!this.m_Dlg) return;
        
        var loadCallback = function(data, txt, rq){
            var formated_html = '';
            this.forEachMarker(function(m) { formated_html += this.parseObjCoords(data, m); }, this);
            if (formated_html.length == 0) { // no markers to display
                formated_html = '<div style="color:gray;margin-top:10px;" align="center">No Markers To Display</div>';
            }
            this.m_Dlg.down('#infoPanel').update(formated_html);
        }
        SeqView.App.simpleAjaxRequest({url: this.getMarkerInfoURL() + '&appname=' + this.m_App.m_AppName, context: this,
        success: loadCallback, error: loadCallback});
    },

    deleteMarker: function(num) {
        var i = this.forEachMarker(function(m) {
            if (m.marker_num == num) {
                m.deleteMarker();
                return false;
            }
        });
        if (i)
            this.updateInfo();
    },

    getMarker: function(mNum) {
        var marker = null;
        this.forEachMarker(function(m) { if (m.marker_num == mNum) marker = m; });
        return marker;
    },

    setMarkerPosition: function(mNum) {
        var marker = null;
        this.forEachMarker(function(tmp) { if (tmp.marker_num == mNum) marker = tmp; });
        if (marker) this.MarkerDlg(null, marker.displayLocalPos(), marker);
    },

    renameMarker: function(num) {
        var m = null;
        this.forEachMarker(function(tmp) {
            if (tmp.marker_num == num) {
                m = tmp;
            }
        });
        if (!m) return; // marker not found

        Ext.MessageBox.prompt('Marker', 'Please enter new marker name:', function(btn, text) {
            if (btn!='ok' || text.length==0) return;
            m.rename(text);
        }, this, false, m.marker_name );
    },

//////////////////////////////////////////////////////////////////////////
// reset:

    reset: function() {
        for(var m in this.m_AllMarkers) {
            if (this.m_AllMarkers[m]) {
                this.m_AllMarkers[m].deleteMarker();
            } // delete markers
        }
        this.updateInfo();
        this.m_AllMarkersRefs = [];
        this.m_MarkerColorIdx = 0;
        this.m_NextMarkerNum = 1;
    },

    updateMarkersPos: function(view) {
        this.forEachMarker(function(m) {
            m.updateMarkerPos(view);
            //console.log('update MarkerPosition');
            //update marker tooltip range text content
            m.m_ToolTip.updateMarkerToolTipContent(this.span, this.seq_pos);
        });
    },

//////////////////////////////////////////////////////////////////////////
// forEachMarker: call fn in scope for each marker, if return of fn is false
//   return number of marker plus one for this function and cancel execution,
//   otherwise return 0.
    forEachMarker: function(fn, scope) {
        var array = this.m_AllMarkersRefs;
        for(var i = 0, len = array.length; i < len; i++){
            var m = array[i];
            if (!m.deleted) {
                if (fn.call(scope || m, m, i, array) === false) { return (i+1); }
            }
        }
        return 0;
    },

    hasMarkers: function() {  return this.m_AllMarkersRefs.length > 0; },
    findMarker: function(m_name) { return this.m_AllMarkers[m_name]; },
    scrollMarkers: function(view, delta) { this.forEachMarker(function(m) { m.scrollPix(view, delta); }); },
    updateMarkersSize: function(view) { this.forEachMarker(function(m) { m.updateMarkerSize(view); }); },

    findMarkerByName: function(name) {
        var array = this.m_AllMarkersRefs;
        for(var i = 0, len = array.length; i < len; i++){
            var m = array[i];
            if (!m.deleted && m.marker_name == name) {
                return m;
            }
        }
        return null;
    },

    findMarkerByPos: function( pos ){
        var array = this.m_AllMarkersRefs;
        for(var i = 0, len = array.length; i < len; i++){
            var m = array[i];
            if( !m.deleted && m.seq_pos == pos ){
                return m;
            }
        }
        return null;
    },

    parseMarkersURL: function(markers, val) {
        if (!val || val.length == 0) return;

        var parts = cbSplit(val, ',');
        while (parts.length) {
            var part = parts.splice(0,1)[0];
            while (parts.length && part.charAt(part.length-1) == '\\') {
                part = part.substr(0, part.length-1) + ',' + parts.splice(0,1)[0];
            }

            var lock = false;
            var lastChar = part.charAt(part.length-1);
            if (lastChar === '!') {
                lock = true;
                part = part.substr(0, part.length-1);
            }

            var info   = part.split('|');
            var coords = info.splice(0,1)[0].split(':');
            var name   = info.splice(0,1)[0];
            while (info.length && name.charAt(name.length-1) == '\\') {
                name = name.substr(0, name.length-1) + '|' + info.splice(0,1)[0];
            }
            name = NCBIGBUtils.unescapeName(name);
            var color  = info.splice(0,1)[0];
            var flags  = info.length > 0 ? parseInt(info.splice(0,1)[0]) : 0;
            if (lock) flags = flags | SeqView.MarkerFlags.UserLock;

            if( isNaN(coords[0]) ) continue;

            var range;
            var marker_pos = NCBIGBUtils.stringToNum(coords[0])-1;
            if( coords.length > 1 ){
                if( isNaN(coords[1]) ) continue;

                range = [marker_pos, NCBIGBUtils.stringToNum(coords[1])-1]
            } else {
                range = [marker_pos];
            }

            var colors = { red: "ff0000", green: "008000", blue: "0000ff", black: "000000", cyan: "00ffff" };
            if (color && color[0] != '#') {
                if (/^[0-9a-fA-F]+$/.test(color)) {
                    color = '#' + color;
                } else if( colors[color] ){
                    color = '#' + colors[color];
                } else {
                    color = '#' + colors['green'];
                }
            }
            markers.push([range, name, flags, color]);
        }
    },

    parseMarkersOldStyle: function(markers, marker_ranges, marker_names) {
        for (var m = 0; m < marker_ranges.length; m++) {
            var val = marker_ranges[m];
            var flags = 0;
            if (val.indexOf('!') != -1) {
                flags = flags | SeqView.MarkerFlags.UserLock;
                val = val.replace(/!/, '');
            }
            var parts = val.split(':');
            var marker_pos = NCBIGBUtils.stringToNum(parts[0])-1;
            var range;
            if (parts.length > 1) {
                range = [marker_pos, NCBIGBUtils.stringToNum(parts[1])-1]
            } else {
                range = [marker_pos];
            }
            markers.push([range, marker_names[m], flags]);
        }
    },

//////////////////////////////////////////////////////////////////////////
// getMarkersData:
// if less_strict == true -- only symbols '|' and ',' will be "escaped"
    getMarkersData: function(less_strict, forPDF) {
        var markers = '';

        this.forEachMarker(function(m) {
            if (m.seq_pos == undefined || isNaN(m.seq_pos) || m.nvMarker) return;
            if (markers.length > 0) markers += ',';

            markers += m.displayPos(forPDF);
            markers += "|" + NCBIGBUtils.escapeName(m.marker_name, less_strict ? /([\|,])/ : false);
            markers += "|" + ((m.color.charAt(0) === '#') ? m.color.slice(1) : m.color);
            if (!forPDF && m.flags) markers += "|" + m.flags;
        });
        if (markers && !forPDF) return '&mk=' + markers;
        return markers;
    },


    getMarkersURL: function() {
        return this.getMarkersData() + (this.m_Dlg && this.m_Dlg.isVisible() ? '&vm=true' : '');
    }
};

Ext.define('SeqView.MarkerToolTip', {
    extend: 'NCBIGBObject.ToolTip',
    width: 180,
    minWidth: 180,
    anchor: 'top',
    initComponent: function() {
        this.updateMarkerToolTipContent(this.selection.span, this.selection.seq_pos);
        this.title = '&nbsp;&nbsp;&nbsp;' + (this.selection.nvMarker ? 'Navigation Marker' : this.selection.marker_name);
        this.callParent(arguments);

        var items = [{xtype: 'tbtext', text: '', itemId: 'range'},
                    {text: (this.selection.nvMarker ? 'Set New Marker at Position' : 'Rename...'), iconCls:'xsv-markers',
                    handler: function() {
                        if (!this.selection.nvMarker) { 
                            this.selection.m_MInfo.renameMarker(this.selection.marker_num);
                            this.pingClick('7-Ren');

                        } else {
                            this.selection.m_MInfo.MarkerDlg(null, this.selection.seq_pos + 1);                            
                            this.pingClick('7-New');
                        }
                    }}];
        if (!this.selection.nvMarker) {
            items.push(
                {text: 'Modify Position/Range/Color', iconCls: 'xsv-zoom_range',
                handler: function() {
                    this.pingClick('7-0');
                    if (!this.isPinned()) this.hide();
                    if (!(this.selection.flags & (SeqView.MarkerFlags.SystemLock | SeqView.MarkerFlags.UserLock))) {
                        this.selection.m_MInfo.setMarkerPosition(this.selection.marker_num);
                    } else 
                        Ext.Msg.alert('Information',
                            (this.selection.flags & SeqView.MarkerFlags.UserLock)
                            ? 'To modify range please first unlock the marker.'
                            : "The object is locked, can't modify range");
            }});
        }
        items.push(
            {text: 'Zoom To Sequence At Marker', iconCls:'xsv-zoom_seq',
            handler: function() {
                this.pingClick('7-1');
                this.selection.zoomSeqMarker(this.gview);
                if (!this.isPinned()) this.hide();
            }},
            {text: 'Marker Details', iconCls: 'xsv-markers',
            handler: function() {
                this.pingClick('7-2');
                if (!this.isPinned()) this.hide();
                this.selection.m_MInfo.showDlg(this.gview);
        }});
       

        var tooltip = this;
        var selection = this.selection;
        var app = this.selection.m_MInfo.m_App;
        if (!app.m_Embedded) {
            var vNum = 0;
            app.forEachView(function(view) {
                if (view.isGraphic() && !(vNum++)) {
                    items.push({
                         text: 'Reveal in Graphic View', iconCls: 'goto_marker', scope: this,
                         handler: function() {
                             this.centerInView(view);
                             this.m_ToolTip.pingClick('7-Graph');
                         }
                    });
                }
            }, selection );
        }
        items.push(
             {text: 'Reveal in Sequence View', iconCls: 'xsv-new_fasta', scope: selection,
             handler: function() {
                 this.showSequence();
                 this.m_ToolTip.pingClick('7-Seq');
             }
        });
        if (!(this.selection.flags & SeqView.MarkerFlags.SystemLock) && !this.selection.nvMarker) {
            items.push(
                {text: 'Lock/Unlock Marker', iconCls: 'xsv-marker_locked',
                handler: function() {
                    this.selection.lockMarker();
                    this.pingClick('7-3');
                }
            });
        }
        items.push(
            {text: 'Remove Marker', iconCls: 'xsv-marker-remove',
            handler: function() {
                if (!this.selection.m_ToolTip.pinned) {
                    this.selection.deleteMarker();
                    this.selection.m_MInfo.updateInfo();
                    this.pingClick('7-4');
                } else {
                    Ext.Msg.alert('Status', 'Can\'t remove marker with pinned tooltip.');
                }
            }
        });

        var scope = this;
        Ext.each(items, function(item, ix){
            item.scope = item.scope || scope; 
            item.style = {margin: '0px', padding: '0px'};
            item.x = 0;
            item.y = 20 * ix;
        });        
        this.toolbar = new Ext.Toolbar({layout: 'absolute', target: this.target, items: items, border: false,
            height: items.length * 20});
       
        this.add(this.toolbar);
        this.addTool([
            { type: 'magnify', scope: this,
                callback: function() {
                    if (!this.isPinned()) this.hide();
                    this.selection.zoomSeqMarker(this.gview);
                }
            }]);
    },

    updateMarkerToolTipContent: function(span, seq_pos) {
        this.adjusted_seq_pos = this.gview.m_App.posToLocal(seq_pos);
        this.adjusted_seq_end = null;
        if (span > 1) {
            this.adjusted_seq_end = this.gview.m_App.posToLocal(seq_pos + span - 1);
            this.rangeText = 'Range: ';
            if (this.gview.m_App.getFlip()) {
                this.rangeText += this.adjusted_seq_end + ' .. ' + this.adjusted_seq_pos;
            } else {
                this.rangeText += this.adjusted_seq_pos + ' .. ' + this.adjusted_seq_end;
            }
            this.rangeText +=  ' (' + span.commify() + ' ' + this.gview.m_App.m_ViewParams.units +')';
        } else
            this.rangeText = 'Position: ' + this.adjusted_seq_pos;
    },

    isPinned: function() { return this.pinned; },

    onTargetOut: function(e){
        if(this.el && !this.disabled && !e.within(this.target.dom, true)){
            this.el.applyStyles('outline-style: none');
        }
        this.callParent(arguments);
    },

    onTargetOver: function(e){
        if (this.isPinned()) {
            if (!this.isVisible()) {
                this.mouseOffset = [0, 0];
                this.show();
            }
            this.el.applyStyles('outline-style: solid;outline-color: red;outline-width:thin');
        }
        else this.callParent(arguments);
//        this.ping({"jsevent":"mouseover", "sv-event":"Marker_ToolTip"});
    },
    onShow: function() {
        if( this.gview.m_ContextMenu && this.gview.m_ContextMenu.isVisible() ){
                this.gview.m_ContextMenu.hide();
        }
        this.callParent(arguments);
        this.rtEl = this.rtEl || this.down('#range');
        this.rtEl.setText('<b><small>' + this.rangeText + '</small></b>');
        this.setWidth(this.rtEl.getWidth() + 10);
        
    },
    ping: function(a) { this.selection.m_MInfo.m_App.ping(a); },
    pingClick: function(a, e) { this.selection.m_MInfo.m_App.pingClick(a, e); }

});
/*  $Id: textview.js 47201 2022-11-02 14:55:11Z borodine $
 * Authors:  Vlad Lebedev, Maxim Didenko
 * File Description: Sequence View Panel
 */

 /********************************************************************/
//////////////////////////////////////////////////////////////////////
// SeqView.TextView
/********************************************************************/


SeqView.TextView = (function() {
    var sm_TextViewRows = 100; // Max of 100 rows in the sequence text panels
    var sm_ViewTmpl = new Ext.Template(
        '<table width="{width}" class="xsv-fasta-table sv-fasta-table" style="border-style:solid;border-width:{border}px;" >',
        '<tbody><tr bordercolor="#ffffff"><td id="left_sequence_{id}" width=50 valign="top" align="right"><br><NOBR><em>{left_nums}</em></NOBR></td>',
        '<td valign="top"><span id="top_sequence_{id}"><em>{top_nums}</em></span><br/><div id="text_sequence_{id}">{sequence}</div></td>',
        '</tr></tbody></table>'
    );

    function constructor(app) {
        var div_id = 'seq-dlg_' + app.m_Idx;
        if (Ext.get(div_id)) { return; }
        this.m_App = app;
        this.m_Idx = 'S' + this.m_App.m_Idx; // m_Idx is needed for Locator object;
        //this.m_View = view;
        this.m_DivId = div_id;
        this.m_Color = '000000';
        this.m_FromSeq = 0;
        this.m_LenSeq = 0;
        this.m_ShowTranslation = true;
        this.m_SequenceCols = Ext.isWindows ? 90 : 100; // 10 rows, 10 bases in each row. 100 bases per line (configurable)
        this.m_SeqStarts = [];
        this.m_Locator = null;
        this.m_PrevXY = [];
        this.m_Flip = app.m_Flip;
        this.m_TextViewSize = this.m_SequenceCols * sm_TextViewRows; // load sequence size (for text views)
        this.m_SelectedFeat = '';
        this.m_HideTrans = 1;

    };

    constructor.getViewTmpl = function() {
        return sm_ViewTmpl;
    };
    constructor.getTextViewRows = function() {
        return sm_TextViewRows;
    };

    return constructor;

})();

////////////////////////////////////////////////////////////////
/// Public methods

SeqView.TextView.prototype = {

    canFlip: function() { return this.m_App.m_ViewParams.acc_type == 'DNA'; },
    isLoading: function() { return this.m_Loading; },
    isPanorama: function() { return false; },
    getURL: function() { return "&seq=" + (this.m_FromSeq+1) + ':' + (this.m_FromSeq+this.m_LenSeq) },
    toSeq: function() { return [this.m_FromSeq, this.m_FromSeq + this.m_LenSeq - 1, this.m_LenSeq]; },
    pingClick: function(a, e) { this.m_App.pingClick(a, e); },

    getTitle: function() {
        var title = "Sequence View";
        if (this.canFlip())
            title += (this.m_Flip ? " (negative" : " (positive") + " strand)";
        return title;
    },

    openDlg: function(obj) {
        if (this.m_SeqDlg) {
            this.m_SeqDlg.expand();
            return;
        }
        var from_seq, len_seq;
        // check if pass a marker position
        if (typeof obj == 'number') {
            from_seq = Math.max(0, obj - 949);
        } else if (Ext.isArray(obj)) {
            // range is passed from the URL
            from_seq = obj[0];
            len_seq = obj[1];
            this.m_Flip = obj[2];
        } else {
            if (!obj.isPanorama) {
                this.m_Flip = obj.orientation == 'reverse';
                this.displaySequenceTextDlg(obj.from, obj.to - obj.from, obj);
                return;
            }
            if (obj.isPanorama()) { // no parent
                from_seq = Math.max(0, this.m_App.m_SeqLength / 2 - this.m_TextViewSize / 2);
            } else {
                this.m_Flip = obj.getFlip();
                var vis_range = obj.toSeq();
                from_seq = Math.max(0, vis_range[0] + vis_range[2] / 2 - this.m_TextViewSize / 2);
            }
            from_seq = Math.round(from_seq/10) * 10;
        }
        len_seq = Math.min(from_seq + this.m_TextViewSize, this.m_App.m_SeqLength) - from_seq;
        this.displaySequenceTextDlg(from_seq, len_seq);
    },

    openDlgHelp:function()
    {
        window.open('/tools/sviewer/sequencepanel/');// + '&inline=true');
    },
    displaySequenceTextDlg: function(from_seq, len_seq) {
        this.m_App.resizeIFrame(500);
        if (typeof from_seq == 'number') {
            this.startImageLoading(from_seq, len_seq);
            if (this.m_App.m_Views[0].isPanorama()) {
                this.m_Locator = new SeqView.Locator(this, this.m_Color, false);
            }
        }
        this.m_SeqDlg = new Ext.Window({
            layout:'fit', title: this.getTitle(),
            app: this.m_App,
            width: 750, height: 450,  minWidth: 480,
            minHeight: 200,
            constrain:false,
            collapsible:true,
            maximizable:true,
            hideAction:'close',
            monitorResize: true,
            cls: 'SeqViewerApp',
            onEsc: function() { this.close(); },
            listeners:{
                'resize': function(window, width, height) {
                    var top_id = 'top_sequence_' + this.m_App.m_Idx;
                    var left_id = 'left_sequence_' + this.m_App.m_Idx;
                    var cur_top = document.getElementById(top_id);
                    if (!cur_top) { return; }
                    var starts_width = Ext.get(left_id).getWidth() + 40;
                    var line_width = cur_top.offsetWidth;
                    var one_width = line_width / this.m_SequenceCols;
                    this.m_SequenceCols = Math.floor((width - starts_width) / one_width / 10) * 10;
                    this.m_TextViewSize = this.m_SequenceCols * SeqView.TextView.getTextViewRows();
                    this.startImageLoading(this.m_FromSeq, Math.min(this.m_TextViewSize, this.m_App.m_SeqLength));
                },
                'maximize': function( window ){
                    window.setWidth(
                        Ext.isIE ? document.documentElement.clientWidth
                        : !Ext.isOpera ? document.body.clientWidth : self.innerWidth
                    );
                },
                'move': function(win, x, y) {
                    var scrollY = Ext.getDoc().getScroll().top;
                    var p = this.m_SeqDlg.el.translatePoints(x, scrollY);
                    if (y < p.top) {
                        this.m_SeqDlg.setPosition(x, p.top);
                    }
                },
                'beforeclose': function() { if (this.m_SeqDlg.maximized) this.m_SeqDlg.restore(); },
                'close': function() {
                    if (this.m_Locator) { this.m_Locator.remove(); }// remove locator form panoram
                    this.m_Locator = null;
                    this.m_App.un('origin_changed', this.onOriginChanged, this );
                    this.m_App.un('strand_changed', this.flipSeqPanelStrand, this );
                    this.m_App.m_TextView = null;
                    delete this.m_SeqDlg;
                },
                scope: this
            },
            id: this.m_DivId,
            cls: 'SeqViewerApp',
            tbar:[
                {text:'Prev Page', scope: this, iconCls:'xsv-prev', handler:this.gotoSequencePrev},
                '-',
                {text:'Next Page', scope: this, iconCls:'xsv-next', handler:this.gotoSequenceNext},
                '->',
                //{text:'Show Gap', iconCls:'xsv-show_gap', handler:SeqApp.sequenceShowGap},'-',
                {xtype: 'combo',
                  itemId: 'tripletCombo',
                  queryMode: 'local',
                  forceSelection: 'true',
                  hidden: true,
                  triggerAction: 'all',
                  store: new Ext.data.ArrayStore({ fields: ['object_id', 'name'] }),
                  displayField: 'name',
                  valueField: 'object_id',
                  width: 250,
                  listeners: { select: function(combo) {
                      this.m_SelectedFeat = combo.getValue();
                      this.startImageLoading(this.m_FromSeq, this.m_LenSeq);
                  }, scope:this }
                },
                ' ',
                {tooltip:'Printer-Friendly Page', iconCls:'xsv-printer', scope:this, handler:function() { SeqView.PrinterFriendlyPage(this.m_App);/*SeqView.App.showPrintPageDlg(this.m_App.m_Idx);*/ }},
                '-',
                {text:'Flip Strands', id:'seq-flip-button', tooltip:'Flip Sequence Strands', iconCls:'xsv-flip-strands',
                    pressed:this.m_Flip, enableToggle:true, hidden: !this.canFlip(),
                    scope: this, handler:function() { this.m_App.setFlip(!this.m_Flip);/*this.flipSeqPanelStrand();*/} },
                this.m_App.m_ViewParams.acc_type != 'DNA' ? {hidden:true} : '-',
                {text: this.translationName(),
                   id: 'seq_trans_btn', //xtype: 'tbsplit',
                   tooltip:'Flip Translation', iconCls:'xsv-show_translation',
                   hidden: this.m_App.m_ViewParams.acc_type != 'DNA',
                   menuAlign:'tr-br?',
                   menu:new Ext.menu.Menu({items:[
                     {text:'Annotated', checked: this.m_HideTrans === 1, group:true, scope: this,
                            handler:function() { this.m_HideTrans = 1; this.flipSeqTranslation(); } },
                     {text:'Conceptual', checked: this.m_HideTrans === 2, group:true, scope: this,
                            handler:function() { this.m_HideTrans = 2; this.flipSeqTranslation(); } },
                     {text:'Both', checked: this.m_HideTrans === 0, group:true, scope: this,
                            handler:function() { this.m_HideTrans = 0; this.flipSeqTranslation(); } },
                     {text:'None', checked: this.m_HideTrans === 3, group:true, scope: this,
                            handler:function() { this.m_HideTrans = 3; this.flipSeqTranslation(); } }
                    ]})
                  },

                  this.m_App.m_ViewParams.acc_type != 'DNA'? {hidden:true} : '-',
                  {text:'Go To Position', iconCls:'xsv-goto_position', scope: this, handler:this.gotoSeqPositionDlg}

            ],
            modal: false,
            hidden: true,
            items:[{ xtype:'panel', autoScroll: true, id:'seq-panel_'+this.m_App.m_Idx, layout:'fit' }],
            buttons: [
                //{xtype: 'tbtext', text: '<a href=\"' + SeqView.webNCBI + 'tools/sviewer/sequencepanel/\" target=\"_blank\">Help</a>'},
                { text:"Help", tooltip:'Help',ariaLabel:'Help', itemId: 'help', dock:'left', scope: this, disabled: false, handler:this.openDlgHelp },
                '->',
            '->',
            {text: 'Close', scope:this, handler: function() { this.m_SeqDlg.close(); }}]
        });
        var app = this.m_App;
        app.on('origin_changed', this.onOriginChanged, this );
        app.on('strand_changed', this.flipSeqPanelStrand, this );
        
        this.m_SeqDlg.on('close', function(){ app.resizeIFrame(); app.m_DialogShown = false; });

        app.m_DialogShown = true;
// This delay is not needed and in fact is buggy - when the callback for startImageLoading is called with the data
// the dialog can be not created yet.
//        if (app.m_iFrame) this.m_SeqDlg.show.defer(500, this.m_SeqDlg); else this.m_SeqDlg.show();
        this.m_SeqDlg.show();
    },

    onOriginChanged: function(app) {
        this.startImageLoading(this.m_FromSeq, this.m_LenSeq);
    },

    startImageLoading: function(from, len) {
        this.m_FromSeq = from;
        this.m_LenSeq = len;
        SeqView.App.simpleAjaxRequest({url: this.getSeqTextURL(), context:this, success: this.seqLoadCallback/*, error: this.seqLoadCallback*/});
    },

    getSeqTextURL: function() {
        var to_seq = this.m_FromSeq + this.m_LenSeq;
        var user_data_key = "";
        if (this.m_App.m_Key) {
            user_data_key = "&key=" + this.m_App.m_Key;
        }
        return this.m_App.m_CGIs.Sequence + '?opts=seq&id=' + this.m_App.GI + '&from=' + this.m_FromSeq + '&to=' + to_seq + '&col=' +
            this.m_SequenceCols + '&translation=' + (this.m_ShowTranslation ? 'true' : 'false') +
            '&reverse=' + (this.m_Flip ? 'true' : 'false') + '&selected=' + this.m_SelectedFeat + user_data_key;
    },

    genTopNumbers: function() {
        var top_nums = "";
        for (var i = 0; i != this.m_SequenceCols / 10; i++) {
            if (this.m_Flip) { top_nums += '<span>0</span>987654321'; }
            else { top_nums += '123456789<span>0</span>'; }
        }
        return top_nums;
    },

    seqLoadCallback: function(data, text) {
        var combo = this.m_SeqDlg.down('#tripletCombo');
        if (!combo) return;
        var from_cgi = data,
            starts = from_cgi.starts;

        if (this.m_App.getOrigin() != 0  ||  this.m_App.getFlip()) {
            var arr,
                new_starts = '',
                re = /([^0-9]*)([0-9]+)(<br>)/ig;
            while ((arr = re.exec(starts)) != null) {
                new_starts += (arr[1] + this.m_App.posToLocal(arr[2]) + arr[3]);
            }
            starts = new_starts;
        }
        var top_nums = this.genTopNumbers();

        var html = SeqView.TextView.getViewTmpl().apply({border:0, width:'95%', left_nums:starts, top_nums:top_nums, id: this.m_App.m_Idx, sequence:from_cgi.sequence});
        Ext.getCmp('seq-panel_'+this.m_App.m_Idx).update(html);
        // Add tooltips
        var idx = 1;
        var ttl_search = from_cgi.app_id? 'ttl-'+from_cgi.app_id+'-' : 'ttl-id';
        var app = this.m_App;
        while (true) {
            var tip_elem = Ext.get(ttl_search + idx);
            if (!tip_elem) { break; }
            var idx_array = tip_elem.dom.getAttribute('rel').split(',');
            var feat_idx = idx_array[idx_array.length-1];
            var sig = "";
            for (j=0; j!=from_cgi.features.length; j++) {
                var feat = from_cgi.features[j];
                if (feat.id == feat_idx) { sig = feat.object_id; break; }
            }
            if (sig.length > 0) {
                var qtip = new Ext.ToolTip({
                    target: tip_elem,
                    dismissDelay: 5000,
                    cfg: { url: app.m_CGIs.ObjInfo+'?objinfo=1&signatures='+sig+'&id='+app.GI,
                        success: function(data) { this.setHtml(data[0].text); }
                    }
                });
                qtip.cfg.context = qtip;
                qtip.on('render', function(){ SeqView.App.simpleAjaxRequest(this.cfg); }, qtip);
            }
            idx++;
        } // while
        // add markers
        this.m_App.forEachMarker( function(m) { this.AddOrUpdateMarker(m); }, this);

        this.flipSeqTranslation();

        if (from_cgi.triplets) {
            combo.getStore().loadData(from_cgi.triplets.triplets);
            this.m_SelectedFeat = from_cgi.triplets.triplets[from_cgi.triplets.selected || 0].object_id;
            combo.setValue(this.m_SelectedFeat);
            combo.show();
            combo.collapse();
        } else {
            combo.hide();
            this.m_SelectedFeat = '';
        }

        Ext.get(this.m_DivId).on({
            'mousedown':    this.onMouseDown,
            'mouseup':      this.onMouseUp,
            'mousemove':    this.onMouseMove,
            'contextmenu':  this.onContextMenu,
            scope: this
        });

        this.m_App.updateLocator(this);
    },


    AddOrUpdateMarker: function(the_marker) {//num, color, seq_pos, name) {
        var seq_elem = Ext.get('text_sequence_' + this.m_App.m_Idx);
        if (seq_elem) {
            var char_height = this.getLineHeight();
            var char_width = this.getCharWidth();
            if (the_marker.span) {
                var marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+the_marker.marker_num+'_b');
                if (!marker) {
                    var tmpl = new Ext.Template('<div orig_id="marker_0_{num}" id="seqtext-marker_{aidx}_{num}_b" data-qtip="start of {name}" class="xsv-text_marker" style="border-color:{color};border-{side}: none; height:{height}; width:{width};"><div data-qtip="Locked" id="seqtext-marker_{aidx}_{num}_lock_b" style="display:{lock};" class="marker_locked_ov"></div> </div>');
                    marker = tmpl.append(seq_elem,
                                    {num:the_marker.marker_num, aidx:this.m_App.m_Idx,
                                     color:the_marker.color, name:the_marker.marker_name,
                                     lock:the_marker.lock?'inline':'none', side:this.m_App.getFlip()?'left':'right'}, true);
                }
                this.setTextMarker(marker, the_marker.seq_pos);
                marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+the_marker.marker_num+'_e');
                if (!marker) {
                    var tmpl = new Ext.Template('<div orig_id="marker_0_{num}" id="seqtext-marker_{aidx}_{num}_e" data-qtip="end of {name}" class="xsv-text_marker" style="border-color:{color};border-{side}: none;height:{height}; width:{width};"><div data-qtip="Locked" id="seqtext-marker_{aidx}_{num}_lock_e" style="display:{lock};" class="marker_locked_ov"></div> </div>');
                    marker = tmpl.append(seq_elem,
                                    {num:the_marker.marker_num, aidx:this.m_App.m_Idx,
                                     color:the_marker.color, name:the_marker.marker_name,
                                     lock:the_marker.lock?'inline':'none', side:this.m_App.getFlip()?'right':'left'}, true);
                }
                this.setTextMarker(marker, the_marker.seq_pos + the_marker.span);
            } else {
                var marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+the_marker.marker_num);
                if (!marker) {
                    var tmpl = new Ext.Template('<div orig_id="marker_0_{num}" id="seqtext-marker_{aidx}_{num}" data-qtip="{name}" class="xsv-text_marker" style="border-color:{color};height:{height}; width:{width};"><div data-qtip="Locked" id="seqtext-marker_{aidx}_{num}_lock" style="display:{lock};" class="marker_locked_ov"></div> </div>');
                    marker = tmpl.append(seq_elem,
                                    {num:the_marker.marker_num, aidx:this.m_App.m_Idx,
                                     color:the_marker.color, name:the_marker.marker_name,
                                     lock:the_marker.lock?'inline':'none'}, true);
                }
                this.setTextMarker(marker, the_marker.seq_pos);
            }
        }
    },


    UpdateLockMarker: function(marker) {
        if (marker.span) {
            var seq_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num+'_lock_b');
            if (seq_marker) {
                seq_marker.setStyle('display', marker.lock ? 'inline' : 'none');
            }
            seq_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num+'_lock_e');
            if (seq_marker) {
                seq_marker.setStyle('display', marker.lock ? 'inline' : 'none');
            }
        } else {
            var seq_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num+'_lock');
            if (seq_marker) {
                seq_marker.setStyle('display', marker.lock ? 'inline' : 'none');
            }
        }
    },


    RemoveMarker: function(marker) {
        if (marker.span) {
            var seq_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num+'_b');
            if (seq_marker)
                seq_marker.remove();
            seq_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num+'_e');
            if (seq_marker)
                seq_marker.remove();
        } else {
            var seq_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num);
            if (seq_marker)
                seq_marker.remove();
        }
    },


    SetSeqPos: function(marker) {
        if (marker.span) {
            var seq_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num);
            if (seq_marker)
                seq_marker.remove();
        } else {
            var seq_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num+'_b');
            if (seq_marker)
                seq_marker.remove();
            seq_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num+'_e');
            if (seq_marker)
                seq_marker.remove();
        }
        this.AddOrUpdateMarker(marker);
    },


    getLineHeight: function() {
        this.txtMetrix = this.txtMetrix || new Ext.util.TextMetrics('top_sequence_' + this.m_App.m_Idx);
        return this.txtMetrix.getHeight("A");
    },
    getCharWidth: function() {
        this.txtMetrix = this.txtMetrix || new Ext.util.TextMetrics('top_sequence_' + this.m_App.m_Idx);
        return this.txtMetrix.getWidth("A");
    },


    setTextMarker: function(marker_elem, seq_pos) {
        if (seq_pos < this.m_FromSeq  ||  seq_pos >= this.m_FromSeq + this.m_LenSeq - 1) {
            marker_elem.setVisible(false);
            return;
        } else marker_elem.setVisible(true);

        var idx,
            left_pos = Ext.get('left_sequence_' + this.m_App.m_Idx).getWidth(),
            brs = Array.from(document.querySelector('#left_sequence_' + this.m_App.m_Idx + ' em').querySelectorAll('br')),
            topline = document.getElementById('top_sequence_' + this.m_App.m_Idx),
            one_width = topline.offsetWidth / this.m_SequenceCols,
            marker_y = - topline.getBoundingClientRect().top - 2;

        // 0-based to 1-based position shift
        seq_pos += 1;

        if (this.m_Flip) {
            idx = brs.findIndex(function(br) { return br.previousSibling.data <= seq_pos; }) - 1;
            if (idx < 0) idx = brs.length - 1;
            left_pos += one_width * (brs[idx].previousSibling.data - seq_pos - 1);
        } else {
            idx = brs.findIndex(function(br) { return br.previousSibling.data > seq_pos; }) - 1;
            if (idx < 0) idx = brs.length - 1;
            left_pos += one_width * (seq_pos - brs[idx].previousSibling.data);
        }
        marker_elem.setLeft(left_pos);
        marker_y += brs[idx].getBoundingClientRect().top;
        marker_elem.setTop(marker_y);
    },


    SeqText2SeqPos: function(new_top, new_left) { // sequence 2 screen (text view)
        var offset = Ext.get('left_sequence_' + this.m_App.m_Idx).getWidth() + 2,
            brs = Array.from(document.querySelector('#left_sequence_' + this.m_App.m_Idx + ' em').querySelectorAll('br')),
            topline = document.getElementById('top_sequence_' + this.m_App.m_Idx),
            one_width = topline.offsetWidth / this.m_SequenceCols,
            start_pos;
        new_top += topline.getBoundingClientRect().top + 5;
        var x = 0,
            idx = brs.findIndex(function(br) { return br.getBoundingClientRect().top > new_top; }) - 1;
        if (idx < 0) idx = brs.length - 1;
        start_pos = parseInt(brs[idx].previousSibling.data) - 1;
        x = Math.min(Math.ceil((new_left - offset) / one_width), this.m_SequenceCols - 1);
        console.log(x, start_pos, x + start_pos, this.m_FromSeq, this.m_TextViewSize, this.m_FromSeq + this.m_TextViewSize);
        return start_pos + (this.m_Flip ? (- x - 1) : x);
    },


    onMouseDown: function(e) {
        if (e.button === 0) {
            var elem = e.getTarget().id;
            if (elem.indexOf('seqtext-marker_') != 0) { return; }

            // check for lock
            var idx = elem.split('_')[2];
            var the_marker = this.m_App.findMarker('marker_0_' + idx);
            if (!the_marker || the_marker.lock) { return; }

            // ok, ready to drag
            this.m_PrevXY = e.getXY();
            this.m_CurMarker = Ext.get(elem);
            this.m_CurMarker.setStyle('cursor', 'move');
        }
    },


    onMouseUp: function(e) {
        if (!this.m_PrevXY || !this.m_CurMarker) { return; }

        var orig = this.m_CurMarker.dom.getAttribute('orig_id');
        var marker = this.m_App.findMarker(orig);

        var seq_pos = this.SeqText2SeqPos(this.m_CurMarker.getTop(true), this.m_CurMarker.getLeft(true));

        if (marker) {
            if (marker.span) {
                var eid = this.m_CurMarker.id;
                if (eid.indexOf('_b') != -1) {
                    var oelem = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num+'_e');
                    if (oelem) {
                        this.setTextMarker(oelem, seq_pos + marker.span -1);
                    }
                    marker.setSeqPos( seq_pos/*-1*/, false );
                } else {
                    var oelem = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num+'_b');
                    var pos = seq_pos - marker.span +1;
                    if (oelem) {
                        this.setTextMarker(oelem, pos);
                    }
                    marker.setSeqPos( pos, false );
                }
            } else {
                marker.setSeqPos( seq_pos/*-1*/, false );
            }
            this.setTextMarker(this.m_CurMarker, seq_pos);
        }

        this.m_PrevXY = null;
        this.m_CurMarker = null;
    },


    onContextMenu: function(e) {
        e.stopEvent();

        var elem = e.getTarget().id;
        //console.log(elem);
        if (elem.indexOf('seqtext-marker') == 0) {
            var menu = new Ext.menu.Menu();
            var idx = elem.split('_')[2];
            var marker = this.m_App.findMarker('marker_0_' + idx);
            if (marker) {
                menu.add(marker.marker_name);
                menu.add('-');
                menu.add({text:'Set To Position...', scope:this, handler:function() {
                    var cur_pos = marker.seq_pos-this.m_App.getOrigin();
                    if (cur_pos >= 0) { cur_pos++; }
                    Ext.MessageBox.prompt('Marker', 'Please enter new sequence position:', function(btn, text) {
                        if (btn!='ok' || text.length == 0) { return; }
                        var seq_pos = this.m_App.convertRelativePosition(text);
                        if (isNaN(seq_pos) || seq_pos < 0 || seq_pos >= this.m_App.m_SeqLength) {
                            Ext.MessageBox.show({title:m.marker_name, msg:'Invalid sequence position.', buttons:Ext.MessageBox.OK, icon:Ext.MessageBox.ERROR});
                            return;
                        }
                        marker.setSeqPos(seq_pos,true);
                    }, this, false, cur_pos ); }
                });
                menu.add('-');
                menu.add({iconCls:'xsv-marker-remove', scope:this, text:'Remove Marker', handler:function() { marker.deleteMarker(); } });
                menu.add({iconCls:'xsv-markers', text:'Marker Details', scope:this, handler:function() { this.m_App.showMarkersDlg( this ); } });

                menu.showAt(e.getXY());
            }
        }
    },


    onMouseMove: function(e) {
        if (!this.m_PrevXY || !this.m_CurMarker) {
            return;
        }

        var delta_x = this.m_PrevXY[0] - e.getXY()[0];
        var delta_y = this.m_PrevXY[1] - e.getXY()[1];
        this.m_PrevXY = e.getXY(); // save new values

        NCBIGBUtils.ClearBrowserSelection();

        var starts_width = Ext.get('left_sequence_' + this.m_App.m_Idx).getWidth() + 2;
        var line_width = document.getElementById('top_sequence_' + this.m_App.m_Idx).offsetWidth;

        var new_left = this.m_CurMarker.getLeft(true) - delta_x;
        var new_top = this.m_CurMarker.getTop(true) - delta_y;

        new_left = Math.min(new_left, starts_width + line_width-7);
        new_left = Math.max(new_left, starts_width);

        new_top = Math.min(new_top, Ext.get('text_sequence_' + this.m_App.m_Idx).getHeight());
        new_top = Math.max(new_top, Ext.get('top_sequence_' + this.m_App.m_Idx).getHeight()+1);

        var seq_pos = this.SeqText2SeqPos(new_top, new_left);

        this.m_CurMarker.setLeft(new_left);
        this.m_CurMarker.setTop(new_top);

        var orig = this.m_CurMarker.dom.getAttribute('orig_id');
        var marker = this.m_App.findMarker(orig);

        if (marker) {
            if (marker.span) {
                var eid = this.m_CurMarker.id;
                if (eid.indexOf('_b') != -1 ) {
                    var oelem = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num+'_e');
                    if (oelem) {
                        this.setTextMarker(oelem, seq_pos + marker.span -1);
                    }
                    marker.setSeqPos( seq_pos/*-1*/, false );
                } else {
                    var oelem = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+marker.marker_num+'_b');
                    var pos = seq_pos - marker.span +1;
                    if (oelem) {
                        this.setTextMarker(oelem, pos);
                    }
                    marker.setSeqPos( pos, false );
                }
            } else {
                marker.setSeqPos( seq_pos/*-1*/, false );
            }
        }
    },


    gotoSequenceNext: function() {
        this.stepSequence(!this.m_Flip);
    },


    gotoSequencePrev: function() {
        this.stepSequence(this.m_Flip);
    },


    stepSequence: function(dir_fwd) {
        var step = dir_fwd ? this.m_TextViewSize : -this.m_TextViewSize;
        var from_seq = Math.min(this.m_FromSeq + step, this.m_App.m_SeqLength - this.m_TextViewSize);
        from_seq = Math.max(0, from_seq);
        this.startImageLoading(from_seq, this.m_LenSeq);
    },


    translationName: function() {
        if (this.m_HideTrans == 1)
            return "Annotated";
        else if (this.m_HideTrans == 2)
            return "Conceptual";
        else if (this.m_HideTrans == 3)
            return "None";
        else
            return "Both";
    },


    flipSeqPanelStrand: function() {
        this.m_Flip = !this.m_Flip;
        this.startImageLoading(this.m_FromSeq, this.m_LenSeq);
        Ext.getCmp('seq-flip-button').toggle(this.m_Flip);
        if (this.m_SeqDlg)
            this.m_SeqDlg.setTitle(this.getTitle());
    },


    flipSeqTranslation: function() {
        var items, i, el;
        if (this.m_HideTrans == 0 || this.m_HideTrans == 3) {
            items = Ext.query('div[class^=xsv-seqtrans_]');
            for(i = 0; i < items.length; i++) {
                el = Ext.get(items[i]);
                el.setStyle('display',(this.m_HideTrans == 0 ? 'block' : 'none'));
                if (this.m_HideTrans == 0) {
                    Ext.util.CSS.updateRule('.xsv-ttu3','border-top-width', '0px');
                    Ext.util.CSS.updateRule('.xsv-ttu7', 'border-top-width', '0px');
                }
            }
        } else {
            Ext.util.CSS.updateRule('.xsv-ttu3', 'border-top-width', '2px');
            Ext.util.CSS.updateRule('.xsv-ttu7', 'border-top-width', '2px');
            items = Ext.query('div[class^=xsv-seqtrans_prot]');
            for(i = 0; i < items.length; i++) {
                el = Ext.get(items[i]);
                el.setStyle('display',(this.m_HideTrans == 1 ? 'block' : 'none'));
            }
            items = Ext.query('div[class^=xsv-seqtrans_trans]');
            for(i = 0; i < items.length; i++) {
                el = Ext.get(items[i]);
                el.setStyle('display',(this.m_HideTrans == 2 ? 'block' : 'none'));
            }
        }
        this.m_App.forEachMarker(function(m) {
            if (m.span) {
                var text_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+m.marker_num+'_b');
                this.setTextMarker(text_marker, m.seq_pos);
                text_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+m.marker_num+'_e');
                this.setTextMarker(text_marker, m.seq_pos + m.span - 1);
            } else {
                var text_marker = Ext.get('seqtext-marker_'+this.m_App.m_Idx+'_'+m.marker_num);
                this.setTextMarker(text_marker, m.seq_pos);
            }
        },this);
        el = Ext.getCmp('seq_trans_btn');
        el.setText(this.translationName());
    },


    gotoSeqPositionDlg: function() {
        Ext.MessageBox.prompt('Go to position',
        'Please enter sequence position<br />(possible&nbsp;value&nbsp;formats&nbsp;are&nbsp;10k,&nbsp;-20,&nbsp;1m):',
            function(btn, text) {
            if (btn!='ok'  || text.length === 0) { return; }
            var position = 0;
            var bad_pos = !NCBIGBUtils.isNumeric(text);
            if (!bad_pos) {
                position = NCBIGBUtils.stringToNum(text);
                position = this.m_App.posToGlobal(position)
                if (position >= this.m_App.m_SeqLength) {
                    bad_pos = true;
                }
            }
            if (bad_pos) {
                Ext.MessageBox.alert('Go to position', 'Invalid position.');
            } else {
                var from_seq = Math.min(position, this.m_App.m_SeqLength - this.m_TextViewSize);
                from_seq = Math.max(0, from_seq);
                this.startImageLoading(from_seq, this.m_LenSeq);
            }
        }, this, false);
    },


    syncToLocator: function() {
        if (!this.m_Locator) {
            return;
        }
        var from_seq = this.m_App.m_Panorama.toSeq( this.m_Locator.getLeft(true) );
        this.startImageLoading(from_seq, this.m_LenSeq);
    }
};

/*  $Id: alignview.js 39699 2017-10-26 20:20:54Z borodine $
 * Authors:  Vlad Lebedev, Maxim Didenko
 * File Description:
 */
 

SeqView.AlignSelection = function(view, area) {
   this.area = area;
   this.m_View = view; // view index
   this.parent_elem = Ext.get(this.m_View.m_DivId);
   
   var tpl = new Ext.Template('<div id="align_selection_id_{idx}" class="over_selection"/>');
   
   this.element = tpl.append(this.parent_elem, {idx:this.m_View.m_Idx}, true);
   this.element.setLeft(area.x); // the output from cgi seems a bit off. That's why this offsets. Or is it browser-specific? :)
   this.element.setTop(area.y); 
   this.element.setHeight(area.h);
   this.element.setWidth(area.w);
   
   this.element.on({
        'click' :   this.m_View.onClick,
        'mousemove' : function(e) {e.stopPropagation();} ,
         scope: this.m_View
   });
   
    this.qtip = new Ext.ToolTip({
        target: this.element, 
        trackMouse:false, 
        autoWidth:true, 
        autoHide:true, 
        html:area.descr, 
        dismissDelay:5000,
        cls: 'SeqViewerApp'
    });
};


SeqView.AlignSelection.prototype = {
    movePix: function(delta) {
        new_left = this.element.getLeft(true) - delta;
        this.element.setLeft(new_left);
    },
    remove: function() { 
        if (this.qtip) {
            this.qtip.destroy();
        }
        this.element.remove(); 
    }
};


/********************************************************************/
//////////////////////////////////////////////////////////////////////
// SeqView.Alignment 
/********************************************************************/

Ext.define('SeqView.Alignment', { 
    extend: 'SeqView.View',
    m_PrevXY: null,
    m_FromSeq: 0,
    m_LenSeq:  0,
    m_ViewPageSize: 500,
    m_Expand: '',
    m_History: [],

    constructor: function(app) { this.callParent(['alignment', app]); },

    isAlignment: function() { return true; },

    createPanel: function() {
        this.m_DivId = 'alignment_id' + this.m_Idx;
        
        var tbar = [];
        
        if (this.m_App.m_Embedded === false) {
            this.m_Spacer = this.m_App.addView({
				border: false, 
				height: this.getSpacerHeight()
			}); // spacer
            tbar.push(this.createChooseColorBtn(), '-');
        }
        var range = this.toSeq();
        var r_range_from = this.m_App.posToLocal(range[0]);
        var r_range_to   = this.m_App.posToLocal(range[1]);
        if( this.m_App.m_Flip ){
            var t = r_range_from;
            r_range_from = r_range_to;
            r_range_to = t;
        }
        
        if( this.m_App.m_Toolbar["history"] == true ){
            tbar.push(
                { iconCls:'back', tooltip: 'Go back', itemId: 'histPrev'+this.m_Idx, scope:this, 
                    handler:function() { 
                        this.stepHistory(); 
                    }
                }
            );          
        }
        
        if( this.m_App.m_Toolbar["name"] == true ){
            var tiptitle = '<b>'
            tiptitle += this.m_App.m_ViewParams['id'];
            tiptitle += ': ' + this.m_App.m_Config.SeqInfo.title;
            tiptitle += '</b><br/><br/>';
            
            tiptitle += r_range_from.commify() + '&nbsp;-&nbsp;' + r_range_to.commify();
             
            tiptitle += '&nbsp;('+ (range[1] - range[0] + 1).commify() + '&nbsp;';
            tiptitle += (this.m_App.m_ViewParams['acc_type']=='protein' ? 'residues' : 'bases');
            tiptitle += '&nbsp;shown';
            
            if( this.m_App.m_ViewParams.acc_type == 'DNA' ){ 
                if( this.m_App.m_Flip ){
                    tiptitle += ',&nbsp;negative strand'; 
                } else { 
                    tiptitle += ',&nbsp;positive strand'; 
                }
            }
            tiptitle += ")";
            
            var title = '<b>' + r_range_from.commify() + ' - ' + r_range_to.commify() + ' (' + range[2].commify() + ' bases shown)</b>';
            
            
            tbar.push({ text: title, itemId: 'tbtitle'+this.m_Idx, tooltip: tiptitle, width: 200, handler: this.openFullView,  scope:this });
        }
        if( this.m_App.m_Toolbar["panning"] == true ){
                tbar.push(
                '-',
                    {iconCls:'pan-left', tooltip: 'Pan left', scope:this, handler:function() { this.gotoPrev(); }},
                    {iconCls:'pan-right', tooltip: 'Pan right', scope:this, handler:function() { this.gotoNext(); }},
                    '-'
                );
            }                       
                            
        if( this.m_App.m_Toolbar["zoom"] == true ){
            tbar.push(
                {text: '&nbsp;&nbsp;-&nbsp;&nbsp;', scale: 'small',tooltip:'Zoom Out', scope:this, handler:function() { this.zoomOut();}}
            );
        
            this.m_TbSlider = new Ext.Slider({
                name: 'zoom', width: 100,  
                minValue: 0, maxValue: 100, value: 0, increment: 5,
                topThumbZIndex: 100,
                tipText: function(thumb){ return String(thumb.value) + '%'; },
                listeners: { 
                    changecomplete: function(el,val) { 
                      var slider_range = 100;
                            var slider_val = 100 - val;
                                      
                            var vis_range = this.toSeq();
                            var zoom_point = vis_range[0] + vis_range[2] / 2;
                                      
                            var max_bases_to_show = this.m_App.m_SeqLength;
                            var rightPix = this.m_View.body.getWidth() - 2;
                            var min_bases_to_show = rightPix * SeqView.MinBpp;
                       
                            var ee = Math.log(max_bases_to_show / min_bases_to_show);
                             
                            var slider_val = Math.min(slider_range, Math.max(0, slider_val));
                             
                            var len = min_bases_to_show * Math.exp(ee * slider_val / slider_range);
                            var from = zoom_point - len / 2;

                            len = Math.min(this.m_App.m_SeqLength, Math.round(len));
                            from = Math.round( Math.max(0, from) );
                            if (from + len > this.m_App.m_SeqLength) {
                                var extra = from + len - this.m_App.m_SeqLength;
                                from -= extra;
                            }
                            this.pushToHistory(range[0], range[1]-range[0],"");
                            this.loadImage(from, len);

                    }.bind(this)
                }
            });
            tbar.push( this.m_TbSlider );   

            tbar.push(
                {text: '&nbsp;&nbsp;+&nbsp;&nbsp;', scale: 'small',tooltip:'Zoom In', scope:this, handler:function() { this.zoomIn();} },
                {iconCls:'xsv-zoom_seq', tooltip:'Zoom To Sequence', scope:this, handler:function() { this.zoomSeq();} }
            );
        }
        tbar.push('->');
        if( this.m_App.m_Toolbar["reload"] == true ){
            tbar.push('-',
                {iconCls:'x-tbar-loading', id:'seq-view-loading-'+this.m_Idx, tooltip:'Reload View', 
                    scope: this, handler:function(){this.refresh(); }, region: 'east'//, style: {"position": "relative", "right":"0"} 
                }
            );
        }
        
        var first_graphic = true;
        for( var i = 0; i < this.m_App.m_Views.length; i++ ){
            var view = this.m_App.m_Views[i];
            if( view && view != this && view.isGraphic() ){
                first_graphic = false;
                break;
            }
        }
        if( first_graphic ){
            if( this.m_App.m_Toolbar["help"] == true ){
            tbar.push(
                {
                    iconCls: 'xsv-question', 
                    tooltip:'Help', 
                    layout: {type:'vbox'},
                    scope: this, 
                    menu:new Ext.menu.Menu({
                        defaultOffsets: [-70,0],
                        scope: this,
                        items:[{
                            text: 'Help',
                            iconCls: 'xsv-question', 
                            scope:this,
                            handler: function() {SeqView.showHelpDlg(); this.pingClick('2-2-10-1');}
                        },{
                            text: 'Link to View', 
                            iconCls: 'xsv-link_to_page',
                            scope: this, 
                            handler: function() {this.m_App.showLinkURLDlg(); this.pingClick('2-2-10-3');}
                        },{
                            text: 'Feedback', 
                            iconCls: 'xsv-feedback',
                            scope: this, 
                            handler: function() {this.m_App.showFeedbackDlg(); this.pingClick('2-2-10-0');}
                        },{
                            text: 'About', 
                            scope: this, 
                            handler: function() {SeqView.showAboutMessage(this.m_Idx); this.pingClick('2-2-10-2');}
                        }]
                    })
                }
            );
            }
         } else {
            tbar.push(
                {iconCls: 'tb-close', tooltip:'Close', scope: this, handler:function() { this.remove(); }}
            );
        }
        
        var tools = [];
        if (this.m_App.m_Embedded === false) {
            tools.push(this.createCloseTBar());
        } else {
            tools.push({id:'help',qtip: 'Help', handler: function(){ SeqView.showHelpDlg(); }})
        }
        
        this.m_View = this.m_App.addView({
            //collapsible:true, title:'Loading...',
            collapsible:false, header: false,
            tools:tools, view:this, tbar:tbar, html:'<div class="alignment_div" id="'+this.m_DivId+'"/>'
        });
        
        
        Ext.get(this.m_DivId).on({
            'mousemove' : this.onMouseMove,
            'contextmenu':this.onContextMenu,
            scope: this
        });
        
        if (this.m_App.m_Embedded === false)
            this.m_Locator = new SeqView.Locator(this, this.m_Color, true);

        this.loadImage(this.m_FromSeq, this.m_LenSeq);
    },
    
    loadImage: function(from,len) {
        this.m_Loading = true; // start loading
        this.m_Width = this.getBodyWidth();
      
        var params = {id: this.m_App.GI, view:'aln', client:'assmviewer', width: this.m_Width, from: from, len:len}
        if (this.m_App.m_AppName && this.m_App.m_AppName.length > 0)
            params.appname = this.m_App.m_AppName;
        var url = this.m_App.m_CGIs.Alignment;
        if (this.m_Expand && this.m_Expand.length>0) params.expand = this.m_Expand;
        this.m_App.addURLParams(params);

        // Add all 'data_key' keys from all 'alignment_track's to the 'key' parameter
        // to feed the data source
        var tracks = this.m_App.getTrackObjects();
        var key = params.key || "";
        for (var i = 0; i < tracks.length; ++i) {
            var track = tracks[i];
            if (track.key === 'alignment_track' && track.data_key) {
                if (key) key += ',';
                key += track.data_key;
            }
        }
        if (key) params.key = key;

      
        Ext.getCmp('seq-view-loading-'+this.m_Idx).disable();
        this.m_App.AjaxRequest({url: url, context: this, data: params,
            success:this.checkJobStatus, error:this.loadFailure});      
    }, 

    checkJobStatus: function(data, text, res) {
        if (!this.m_Loading) return;

        var from_cgi = SeqView.decode(data);
        if (from_cgi.job_status) {
            if (from_cgi.job_status == 'failed') {
                this.loadFailure(null, from_cgi.error_message);
            } else if(from_cgi.job_status == 'canceled') {
                this.loadFailure(null, 'Job canceled');
            } else {
                var url = this.m_App.m_CGIs.Alignment + '?job_key=' + from_cgi.job_id
                Ext.defer(this.m_App.AjaxRequest, 2000,this,[{url:url, context: this,
                        success: this.checkJobStatus, error: this.loadFailure}]);
            }
        } else {
            Ext.getCmp('seq-view-loading-'+this.m_Idx).enable();
            if (from_cgi.error) {
                this.loadFailure(null, from_cgi.error);
            } else if (from_cgi.success === false) {
                this.loadFailure(null, from_cgi.msg);
            } else {
                // If the img_url begins with ? it contains only parameters for ncfetch, so prepend ncfetch URL
                // This is a way to provide reliable URL resolution for embedding. SV-1760
                if (from_cgi.img_url && from_cgi.img_url.charAt(0) == '?') {
                    from_cgi.img_url = this.m_App.m_CGIs.NetCache + from_cgi.img_url;
                }
                this.m_FromCgi = from_cgi;
                this.m_FromSeq =      this.m_FromCgi.from; // got from cgi 0-based
                this.m_LenSeq =       this.m_FromCgi.len;
                this.m_AlgnLen =      this.m_FromCgi.seq_length;
                this.m_ViewPageSize = this.m_LenSeq;
                this.m_Width =        this.m_FromCgi.img_width;
                this.m_Height =       this.m_FromCgi.img_height;
                this.m_Loading =      false; // loaded
    
                var the_div = Ext.get(this.m_DivId);
                var img_el = the_div.first('img');
                if (img_el) {
                    var d = Ext.getDom(img_el);
                    d.src = from_cgi.img_url;
                } else {
                    var tpl = new Ext.Template('<img class="sv-drag sv-highlight sv-dblclick" src="{img_url}">'),
                    img_el = tpl.append(the_div, from_cgi,true);
                }
                the_div.setStyle('height', this.m_Height + 'px');
                this.m_View.updateLayout();
                if( this.m_TbSlider ){  
                    var slider_range = 100;

                    var range = this.toSeq();
                    this.pushToHistory(range[0], range[1]-range[0],"");
                    var vis_len = range[2];
                              
                    var max_bases_to_show = this.m_App.m_SeqLength;
                    var rightPix = this.m_View.body.getWidth() - 2;
                    var min_bases_to_show = rightPix * SeqView.MinBpp;
                    //var min_bases_to_show = this.m_App.getPanoramaWidth() / 10;
                              
                    var ee = Math.log(max_bases_to_show / min_bases_to_show);
                     
                    var slider_val = slider_range * Math.log(vis_len / min_bases_to_show) / ee;

                    this.m_TbSlider.setValue( 100 - slider_val );
                }
                this.m_App.notifyViewLoaded(this);
            }            
        }
    },
    loadFailure: function(data, res) {
        Ext.MessageBox.show({title:'Image loading error', msg:res, buttons:Ext.MessageBox.OK, icon:Ext.MessageBox.INFO});
        var the_div = Ext.get(this.m_DivId);
        the_div.setStyle('background-image', 'none');
        this.m_View.setTitle("Not available");
    },

    onMouseDown: function(e) {},
    onMouseUp: function(e) {},
    onMouseMove: function(e) {
        var area = this.hitTest(e.getXY());
        
        var elem = Ext.get( this.m_DivId );
        
        if (area && !this.m_Selection) {
            elem.setStyle('cursor', 'pointer');
            this.m_Selection = new SeqView.AlignSelection(this, area);
        } else if (this.m_Selection) {
            elem.setStyle('cursor', 'default');
            this.m_Selection.remove();
            this.m_Selection = null;
        }
    },

    onClick: function(e)
    {
        var area = this.hitTest(e.getXY());
        if (area) { // open alignment
            var cur_name = area.name;
            var cur_expand = this.m_Expand;
            var e_array = cur_expand.length > 0 ? cur_expand.split(',') : [];
            var expanded = false;

            for (var i = 0;  i != e_array.length;  i++) {
                if (e_array[i] == cur_name) {
                    expanded = true;
                    delete e_array[i];
                    break;
                }
            }
      
            if (!expanded) e_array.push(cur_name);
            this.m_Expand = e_array.join(',');      
            this.loadImage(this.m_FromSeq, this.m_LenSeq);
        }
    },

    onMouseOut: function(e) {
    },

    onContextMenu: function(e) {
        e.stopEvent();
        menu = new Ext.menu.Menu();
             
        menu.add({iconCls:'xsv-zoom_plus', text:'Zoom In', scope:this, handler:function() { this.zoomIn();} });
        menu.add({iconCls:'xsv-zoom_minus', text:'Zoom Out', scope:this, handler:function() { this.zoomOut();} });
        menu.add('-');
        menu.add({iconCls:'xsv-zoom_seq', text:'Zoom To Sequence', scope:this, handler:function() { this.zoomSeq();} });
        menu.add('-');
        menu.add({iconCls:'expand_all', text:'Expand All', scope:this, handler:function() { this.expandAll();} });
        menu.add({iconCls:'collapse_all', text:'Collapse All', scope:this, handler:function() { this.collapseAll();} });
        menu.add('-');
        menu.add({iconCls:'prev', text:'Prev Page', scope:this, handler:function() { this.gotoPrev();} });
        menu.add({iconCls:'next', text:'Next Page', scope:this, handler:function() { this.gotoNext();} });
        menu.showAt(e.getXY());
    },

    toSeq: function() { 
         return [this.m_FromSeq, this.m_FromSeq + this.m_LenSeq - 1, this.m_LenSeq];
    },
    
    updateTitle: function() {
        var range = this.toSeq();
        var r_range_from = range[0] + 1;
        var r_range_to   = range[1] + 1;

        var title = this.m_App.m_Embedded !== false ? "<a href='#' class='hidden_for_print' id='new_view_link_al"+this.m_App.m_Idx+"' style='float:right;padding-right:7px;'>Open Full View</a>" : 'Alignment View: ';
        title += r_range_from.commify() + ' - ' + r_range_to.commify() + ' (' + range[2].commify() + ' bases shown)';
        this.m_View.setTitle(title);

        var href = Ext.get('new_view_link_al'+this.m_App.m_Idx); 
        if (href) { 
            href.on({ 'click' : this.openFullView,  scope:this });
        }

        this.m_App.updateLocator(this);
    },

    hitTest: function(page_xy) {
        var elem_xy = Ext.get( this.m_DivId ).getXY();
        var xx = page_xy[0] - elem_xy[0];// - config['scroll_pix'];
        var yy = page_xy[1] - elem_xy[1];// - config['top_offset'];
        for (var a = 0; a < this.m_FromCgi.checkboxes.length; a++) {
            var area = this.m_FromCgi.checkboxes[a];
            var the_x = area.x;
            var the_y = area.y;
            var the_w = area.w;
            var the_h = area.h;
            if (xx >= the_x-1 && xx <= the_x+the_w && yy >= the_y-1 && yy <= the_y+the_h) return area; // add extra pixel
        }
        return null;
    },
    
    zoomIn: function() {
    	 var new_len  = this.m_LenSeq / 2;
        var new_from = this.m_FromSeq + new_len / 2;
        this.loadImage(new_from, new_len);
    },

    zoomSeq: function( center_seq_pos ){        
        
        if( !center_seq_pos ){
            center_seq_pos = Math.floor( this.m_FromSeq + this.m_LenSeq /2 );
        }

        var new_len  = Math.floor( this.getBodyWidth() * SeqView.MinBpp ); 
        var new_from = Math.floor( center_seq_pos - new_len /2 );
        if( new_from < 0 ){
            new_from = 0;
        }
        if( new_from + new_len > this.m_App.m_SeqLength ){
            new_len = this.m_App.m_SeqLength - new_from;
        }

        this.loadImage( new_from, new_len );
    },


    zoomOut: function() {
        var new_len  = this.m_LenSeq * 2;
        var new_from = this.m_FromSeq - this.m_LenSeq / 2;
   
        new_from = Math.max(0, new_from); // not to exteed the sequence range
        if (new_from + new_len > this.m_AlgnLen) new_len = this.m_AlgnLen - new_from;

        this.loadImage(new_from, new_len);
    },
    

    gotoPrev: function() {
        var new_from = Math.max(0, this.m_FromSeq - this.m_ViewPageSize);
        this.loadImage(new_from, this.m_LenSeq);
    },
    

    gotoNext: function() {
        var new_from = Math.min(this.m_FromSeq  + this.m_ViewPageSize - 1, this.m_AlgnLen - this.m_ViewPageSize);
        new_from = Math.max(0, new_from);
        this.loadImage(new_from, this.m_LenSeq);
    },  

    expandAll: function() {
        var e_array = [];
        for (var a = 0; a < this.m_FromCgi.checkboxes.length; a++) 
            e_array.push( this.m_FromCgi.checkboxes[a].name );
   
        this.m_Expand = e_array.join(',');   
        this.loadImage(this.m_FromSeq, this.m_LenSeq);
    },

    collapseAll: function() {
        this.m_Expand = '';
        this.loadImage(this.m_FromSeq, this.m_LenSeq);
    },
    
    refresh: function(options) {
        if (this.m_VisLenSeq == 0)
            return;
        this.loadImage(this.m_FromSeq,this.m_LenSeq);
    },


    syncToLocator: function() {
         if (!this.m_Locator) {
            return;
        }            
        var vis_from, vis_to;
        if (this.m_Locator.m_ResizeRight) { vis_from = this.m_FromSeq; }
        else { vis_from = this.m_App.m_Panorama.toSeq( this.m_Locator.getLeft(true) ); }

        if (this.m_Locator.m_Scroll) { vis_to = vis_from + this.m_LenSeq-1; }// keep length 
        else { vis_to = this.m_App.m_Panorama.toSeq( this.m_Locator.getLeft(true) + this.m_Locator.getWidth()+3 ); }
        this.loadImage(vis_from, vis_to-vis_from+1);
    },
    
    checkLocatorWidth: function(width) { return true; },
    
    openFullView: function() { this.m_App.getLinkToThisPageURL('portal'); },
    
    pushToHistory: function( from, length, title ){
        if (this.m_HistIx < 0) this.m_HistIx = 0;
        var current = { from: from, len: length, title: title };

        if (this.m_HistIx < this.m_History.length) {
            if (this.mf_HistoricalBack
                || (Math.abs(current.from - this.m_History[this.m_HistIx].from) <= 1
                    && current.len == this.m_History[this.m_HistIx].len)) {
                this.mf_HistoricalBack = false;
                return;
            }

            if( this.m_HistIx > 0 ){
                this.m_History = this.m_History.slice( this.m_HistIx );
                this.m_HistIx = 0;
            }
            
            this.m_History.unshift( current );
            
        } else {
            this.m_History.push( current );
            this.m_HistIx = this.m_History.length -1;
        }
        this.mf_HistoricalBack = false;
        
        var MAX_HISTORY_LENGTH = 10;
        var mhl = MAX_HISTORY_LENGTH;
        
        if( this.m_History.length > mhl ){
            if( this.m_HistIx < mhl ){
                this.m_History = this.m_History.slice( 0, mhl );
            } else {
                this.m_History = this.m_History.slice( this.m_HistIx -mhl +1, mhl );
                this.m_HistIx = this.m_History.length -1;
            }
        }
        this.updateHistoryButtons();
    },


    stepHistory: function( forward ){
        if( forward ){
            if( this.m_HistIx > 0 ){
                this.m_HistIx--;
                this.mf_HistoricalBack = true;
                this.loadImage( this.m_History[this.m_HistIx].from, this.m_History[this.m_HistIx].len ); 
            }
        } else {
            if( this.m_HistIx < this.m_History.length-1 ){
                this.m_HistIx++;
                this.mf_HistoricalBack = true;
                this.loadImage( this.m_History[this.m_HistIx].from, this.m_History[this.m_HistIx].len ); 
            }
        }
        this.updateHistoryButtons();
    },  


    updateHistoryButtons: function() {
        var tbar = this.m_View.getDockedComponent(0);
        if( tbar ){
            var btn_prev = tbar.getComponent( 'histPrev'+this.m_Idx );
            if( btn_prev ){
                var tooltip = '';
                if( this.m_HistIx >= this.m_History.length-1 ){
                    btn_prev.disable();
                    tooltip = 'Back';

                } else {
                    btn_prev.enable();
                    
                    if( this.m_HistIx+1 < this.m_History.length ){
                        tooltip += 'Back to ' + this.m_History[this.m_HistIx+1].title;
                    }
                }
                btn_prev.setTooltip( tooltip );
            }
            
            var btn_next = tbar.getComponent( 'histNext' );
            if( btn_next ){
                var tooltip = '';
                if( this.m_HistIx <= 0 ){
                    btn_next.disable();
                    tooltip = 'Forward';

                } else {
                    btn_next.enable();
                    
                    for( var j = this.m_HistIx-1; j >= 0; j-- ){
                        tooltip += this.m_History[j].title + '<br/>';
                    }
                }
                btn_next.setTooltip( tooltip );
            }
        }
    }

});

/* $Id: graphview.js 47494 2023-05-15 14:18:52Z borodine $
 * Authors:  Vlad Lebedev, Maxim Didenko, Victor Joukov
 * File Description:
 */

SeqView.ChunkWidth = 4000;
SeqView.MinBpp = 1/24;
SeqView.PAN_RIGHT = 1;
SeqView.PAN_LEFT = -1;

/********************************************************************/
//////////////////////////////////////////////////////////////////////
// SeqView.Graphic
/********************************************************************/

Ext.define('SeqView.Graphic', {
    extend: 'SeqView.View',
    m_trackTag: Ext.isFirefox ? 'img' : 'canvas',
    m_bottomRulerHeight: 0,
    m_FromSeq: 0,
    m_LenSeq:  0,
    m_VisFromSeq: 0,
    m_VisLenSeq: 0,
    m_ScrollPix: 0,
    m_Flip: false,
    m_Selection: null,
    m_RangeSelectionSet: [],
    m_SelectedSig: null,
    m_SelectedRangeSet: [],
    m_TbSlider: null,
    m_showAllGenes: false,


    constructor: function(app) {
        this.clear();
        SeqView.Graphic.superclass.constructor.apply(this, ['graphical',app]);
        this.m_LenSeq = this.m_App.m_SeqLength;
        this.m_Flip = app.getFlip();
        if (this.m_App.m_hiddenOptions.indexOf('use_canvas') >= 0) this.m_trackTag = 'canvas';
        if (this.m_App.m_hiddenOptions.indexOf('use_img') >= 0) this.m_trackTag = 'img';
    },

    isGraphic: function() { return true; },
    canFlip: function() { return this.m_App.m_ViewParams.acc_type == 'DNA'; },
    flipStrand: function() { this.setFlip(!this.getFlip()); },
    getFlip: function() { return this.m_Flip; },
    getHeight: function() { return this.bodyDiv.getHeight(); },
    getMinLength: function(){ return Math.floor(this.getBodyWidth() * SeqView.MinBpp); },
    getMostRightPix: function() { return this.getBodyWidth(); },
    getScreenHeight: function() { return this.m_View.body.getHeight() - (this.m_App.m_Embedded ?  0 : 2); },

    setFlip: function(flip) {
        if (flip != this.m_Flip) {
            this.setFlipNoReload(flip);
            this.refresh();
        }
    },

    setFlipLocal: function( flip ){
        if( flip != this.m_Flip ){
            this.m_Flip = flip;
            var button = Ext.getCmp( 'flip-button-' + this.m_Idx );
            if( button ) button.toggle( flip );
            this.refresh();
            this.pushToHistory();
        }
    },

    setFlipNoReload: function(flip,no_broadcast) {
        this.m_Flip = flip;
        var button = Ext.getCmp('flip-button-' + this.m_Idx);
        if (button) button.toggle(flip);
        if (typeof no_broadcast == "undefined" || !no_broadcast) {
            // Propagate to the peer views
            var this_view = this;
            this.m_App.setFlip(flip, this_view);
        }
    },

    seq2Pix: function(seq_pos) {
        var pos = 0;
        if (this.getFlip()) {
            pos = (this.m_LenSeq - (seq_pos - this.m_FromSeq) - 1 + 0.5) * this.m_bpPix;
        } else {
            pos = (seq_pos - this.m_FromSeq + 0.5) * this.m_bpPix;
        }
        return Math.round(pos);
    },

    seq2PixScrolled: function(seq_pos) { return this.seq2Pix(seq_pos) + this.m_ScrollPix; },

    pix2Seq: function(pix_pos, adj) { // screen 2 sequence (gview)
        if (adj === undefined) adj = 0;
        if (this.getFlip()) pix_pos = this.m_Width - pix_pos;
        return  this.m_FromSeq + Math.floor(pix_pos / this.m_bpPix + adj);

    },

    getTrackByArea: function(area) {
        var sig = area.signature.split(';')[0],
            trk = this.m_graphTracks.find(function(t) { return t.signature == sig || t.id == area.parent_id; } );
        return trk ? this.m_App.m_Config.TrackConfig[trk.cfgIdx] : {};
    },

    toSeq: function() { // screen 2 sequence range (more precise calculation required!!!)
        // If we positioned exacly as requested and did not move,
        // use exact cached position
        if (this.m_ReqFrom !== undefined && this.m_ReqLen !== undefined)
            return [this.m_ReqFrom, this.m_ReqFrom + this.m_ReqLen - 1, this.m_ReqLen];

        var ppb = this.m_bpPix;
        var view_seq_width = Math.floor(this.getBodyWidth()/ppb);

        var from, to;
        if (this.getFlip()){
            var left = this.m_FromSeq + this.m_LenSeq - 1;
            to = left - Math.round(-this.m_ScrollPix/ppb);
            from = to - view_seq_width + 1;
        } else {
            from = this.m_FromSeq + Math.round(-this.m_ScrollPix/ppb);
            to = from + view_seq_width -1;
        }

        if (to >= this.m_App.m_SeqLength) to = this.m_App.m_SeqLength - 1;

        return [from, to, to - from + 1];
    },


    createPanel: function() {
        var modeButtHandler = function(b) {
            if (b.iconCls.indexOf('xsv-mode') == 0) {
                this.m_slimMode = b.pressed;
                b.setIconCls(b.pressed ? 'xsv-mode_orig' : 'xsv-mode_slim');
                b.setTooltip('Switch to ' + (b.pressed ? 'normal' : 'slim') + ' mode');
                this.pingClick('2-2-11');
            } else {
                this.setGeneMode(b.pressed);
                this.pingClick('2-2-GeneMode');
            }
            this.refresh();
        };
        var app = this.m_App;

        var slimModeButton = {iconCls: this.m_slimMode ? 'xsv-mode_orig' : 'xsv-mode_slim', scope: this,
            tooltip: 'Switch to ' + (this.m_slimMode ? 'normal' : 'slim') + ' mode',
            pressed: this.m_slimMode, enableToggle: true, handler: modeButtHandler};
        var geneModeButton = {iconCls: 'xsv-gene_mode', itemId: 'btnGeneMode', scope: this,
            pressed: false, enableToggle: true, handler: modeButtHandler};
        var reloadButton = {iconCls:'x-tbar-loading', disabledCls: 'xsv-search-loading', scope: this,
            hidden: !(app.m_Toolbar["reload"] == true),
            itemId: 'btnReload', tooltip: 'Reload View', ariaLabel: 'Reload View',
            handler:function(b, e){ this['refresh' + (e.ctrlKey ? 'CTRL' : '')]({buttonAction: true}); this.pingClick('2-2-9');}};
//            handler:function(){ this.refresh({buttonAction: true}); this.pingClick('2-2-9');}};
        var seqInfo = app.m_ViewParams;
        this.m_DivId = 'graphical_id' + this.m_Idx;
        this.ping({"embedded": app.m_Embedded});
        if (app.m_Embedded && app.m_Embedded == 'minimal') {
            this.m_View = app.addView({
                layout:'hbox',
                header: false,
                border: false,
                align: 'stretchmax',
                items: [{
                    xtype: 'buttongroup',
                    cls: 'hidden_for_print',
                    columns: 1,
                    width: 24,
                    frame: false,
                    items: [
                        {iconCls: 'xsv-full_view', tooltip: 'Switch to full standalone view', scope:this, 
                         handler: function() {this.openFullView('full'); this.pingClick('2-2-2-2');}},
                        {iconCls:'xsv-flip-strands', id:'flip-button-'+this.m_Idx, tooltip:'Flip Sequence Strands',
                         pressed: this.getFlip(), enableToggle:true, hidden: !this.canFlip(), scope: this,
                         handler:function() { this.flipStrand(); this.pingClick('2-2-7-2');;}},
                        {iconCls: 'xsv-zoom_plus', tooltip:'Zoom In', scope:this,
                         handler:function() {this.zoomIn(); this.pingClick('2-2-5-1');}},
                        {iconCls: 'xsv-zoom_minus', tooltip:'Zoom Out', scope:this,
                         handler:function() {this.zoomOut(); this.pingClick('2-2-5-1');}},
                        {iconCls: 'xsv-zoom_seq', tooltip:'Zoom To Sequence', scope:this,
                         handler:function() {this.zoomSeq(); this.pingClick('2-2-6');}},
                        {iconCls:'xsv-goto_position', tooltip:'Go To Position/Range', scope: this,
                         handler:function() {this.gotoPositionDlg(); this.pingClick('2-2-3');}},
                        geneModeButton,
                        slimModeButton,
                        reloadButton,
                        {iconCls: 'xsv-config', tooltip:'Configure tracks', ariaLabel :'Configure tracks',scope: this, hidden: (app.m_NoConfDlg === true),
                         handler:function() { app.showTracksConfigDlg('2-2-8-0'); }},
                        {iconCls: 'xsv-question', tooltip:'Help',ariaLabel:'Help', scope: this,
                         handler:function() {SeqView.showHelpDlg(); this.pingClick('2-2-10');}}
                    ]
                },{
                    header: false,
                    border: false,
                    view:this,
                    flex: 1,
                    autoHeight: true,
                    layout: 'fit',
                    html: {
                        tag: 'div', id: this.m_DivId,
                        cls: 'graphical_div sv-drag sv-highlight sv-dblclick'
                    }
                }]
            });
        } else {
            var tbar = [];
            if (app.m_Embedded === false) {
                this.m_Spacer = this.m_Spacer || app.addView({
                    border: false,
                    height: this.getSpacerHeight()
                });
                if (!this.m_Color || this.m_Color.length == 0) {
                    var cp = new Ext.ColorPalette();
                    this.m_Color = cp.colors[Math.floor( Math.random() * cp.colors.length)];
                }
            }

            var r_range_from = app.posToLocal(this.m_FromSeq);
            var r_range_to   = app.posToLocal(this.m_FromSeq + this.m_LenSeq - 1);
            if (this.getFlip()) {
                var t = r_range_from;
                r_range_from = r_range_to;
                r_range_to = t;
            }

            if (app.m_Toolbar["history"] == true) {
                tbar.push(
                    {iconCls: 'back', ariaLabel: 'Go back', tooltip: 'Go back', itemId: 'histPrev', scope:this,
                     handler:function() {this.stepHistory(false, {from_ui: true}); this.pingClick('2-2-1');}}
                );
            }

            if (app.m_Toolbar["name"] == true) {
                var tiptitle = '<b>' + seqInfo.id + ': ' + seqInfo.title + '</b><br><br>' + r_range_from.commify()
                + '&nbsp;-&nbsp;' + r_range_to.commify() + '&nbsp;(' + this.m_LenSeq.commify() + '&nbsp;'
                + seqInfo.units + '&nbsp;shown';

                if (this.canFlip())
                    tiptitle += ',&nbsp;' + (this.getFlip() ? 'negative' : 'positive') + 'strand';
                tiptitle += ")";

                var menu_items = [
                    {text: 'View GenBank record', tooltip: 'Go to GenBank flat file', scope: this,
                     handler: function() {this.openFullView('portalGB'); this.pingClick('2-2-2-1');}},
                    {text: 'View in GDV', tooltip: 'Go to Genome Data Viewer', scope: this, hidden: true, itemId: 'GDV',
                     handler: function() {this.openFullView('GDV'); this.pingClick('2-2-2-GDV');}},
                    {text:'Full View', tooltip: 'Switch to full standalone view', scope: this, iconCls: 'xsv-full_view',
                     handler: function() {this.openFullView('full'); this.pingClick('2-2-2-2');}}
                ];
                if (seqInfo.assm_info) {
                    menu_items.push({itemId: 'assmInfo', text: 'Assembly ' + seqInfo.assm_info.name, scope: this,
                        tooltip: 'Assembly context information',
                        ariaLabel: 'Assembly context information',
                        icon: seqInfo.icon,
                        handler: function() { this.showAssmInfo(); this.pingClick('2-2-2-3'); }
                    });
                }
                if (seqInfo.placements_status == 4)
                    menu_items.push('-', {text: 'Loading...', iconCls:'xsv-search-loading', disabled: true, itemId: 'molPlacement'});
                else
                    this.checkMolPlacement = Ext.emptyFn;
                

                tbar.push({
                    text: '<b style="font-size:x-small;">' + seqInfo.id + '</b>',
                    tooltip: tiptitle,
                    ariaLabel: seqInfo.id,
                    itemId: 'btnTitle',
                    iconCls: 'xsv-seq-logo',
                    menu: new Ext.menu.Menu({ items: menu_items }),
                    listeners: { menushow: function(btn, menu)
                        { if (this.m_App.m_ViewParams.mappedGDV) menu.down('#GDV').setHidden(false);
                    }, scope: this},
                    scope: this
                });
            }
            if (app.m_Toolbar["search"] == true) {
                var tt = '<b>Go to position/range:</b><br>'
                        + 'Range formats are 10k-20k, -20--10, -10k:-5, 5 to 515, -1m..1m <br><br>'
                        + '<b>Search:</b><br>'
                        + 'Feature name, component name, HGVS, SNP rs id, track info, or sequence (nucleotide regexp with IUPAC equivalents, PROSITE patterns)';
                var buttonSGT = new Ext.Button({
                    text : 'Find:',
                    ariaLabel: 'Find:',
                    itemId: 'btnFind',
                    scope: this,
                    handler: function() {this.gotoAndSearch(); this.pingClick('2-2-3');},
                    tooltip: tt
                });
                var gotoBox = Ext.create('Ext.form.field.ComboBox', {
                    itemId: 'gotoBox', width:150,
                    ariaLabel: 'find',
                    hidden: true, queryModel: 'local', store: app.searchPatternData,
                    valueField: 'pattern', displayField: 'pattern', typeAhead: true,
                    listeners: {
                        specialkey: function(f,e){
                            if (e.getKey() == e.ENTER){
                               this.gotoAndSearch();
                               this.pingClick('2-2-3');
                               e.stopEvent();//to stop propagation affecting other non-SV text fields in the embedded sviewer
                            }
                        }.bind(this)
                    }
                });

                tbar.push(seqInfo.acc_type != 'DNA'? {hidden:true} : '-', buttonSGT, gotoBox);
                //smaller footprint of the above "Find on Sequence" button and its text field for cases when browser's width is
                //less than 900 px
                var buttonSearch = new Ext.Button({
                    itemId: 'btnSearch',
                    hidden: true,
                    tooltip: tt,
                    iconCls: 'xsv-search-button',
                    handler: function() { app.showSearchParamsDlg(this); this.pingClick('2-2-3');},
                    scope: this
                });
                tbar.push(buttonSearch,'-');
            }

            if( app.m_Toolbar["panning"] == true ){
                tbar.push(
                    {iconCls:'pan-left', tooltip: 'Pan left', ariaLabel: 'Pan left', scope:this,
                     handler:function() {this.scrollViewTo(this.m_ScrollPix + this.getBodyWidth()/2); this.pingClick('2-2-4');}},
                    {iconCls:'pan-right', tooltip: 'Pan right',ariaLabel: 'Pan right', scope:this,
                     handler:function() {this.scrollViewTo(this.m_ScrollPix - this.getBodyWidth()/2); this.pingClick('2-2-4');}},
                    '-'
                );
            }

            if( app.m_Toolbar["zoom"] == true ){
                tbar.push(
                    {iconCls:'xsv-zoom_minus', tooltip:'Zoom Out',ariaLabel:'Zoom Out', scope: this,
                     handler: function() {this.zoomOut(); this.pingClick('2-2-5-1');}}
                );

                this.m_TbSlider = new Ext.Slider({
                    name: 'zoom', width: 100,
                    minValue: 0, maxValue: 100, value: 0, increment: 5,
                    topThumbZIndex: 100,
                    tipText: function(thumb){ return String(thumb.value) + '%'; },
                    listeners: {
                        changecomplete: function(el,val) {
                            this.pingClick('2-2-5');
                            var slider_range = 100;
                            var slider_val = 100 - val;

                            var vis_range = this.toSeq();
                            var zoom_point = vis_range[0] + vis_range[2] / 2;

                            var max_bases_to_show = app.m_SeqLength;
                            var min_bases_to_show = this.getMostRightPix() * SeqView.MinBpp;

                            var ee = Math.log(max_bases_to_show / min_bases_to_show);

                            var slider_val = Math.min(slider_range, Math.max(0, slider_val));

                            var len = min_bases_to_show * Math.exp(ee * slider_val / slider_range);
                            var from = zoom_point - len / 2;

                            len = Math.min(app.m_SeqLength, Math.round(len));
                            from = Math.round( Math.max(0, from) );
                            if (from + len > app.m_SeqLength) {
                                var extra = from + len - app.m_SeqLength;
                                from -= extra;
                            }
                            this.startImageLoading(from, len, {from_ui: true});

                        }.bind(this)
                    }
                });
                tbar.push( this.m_TbSlider );

                tbar.push(
                    {iconCls:'xsv-zoom_plus', tooltip: 'Zoom In',ariaLabel: 'Zoom In', scope: this,
                     handler:function() {this.zoomIn(); this.pingClick('2-2-5-1');}},
                    {iconCls: 'xsv-zoom_seq', tooltip: 'Zoom To Sequence',ariaLabel: 'Zoom To Sequence', scope: this,
                     handler:function() {this.zoomSeq(); this.pingClick('2-2-6');}}
                );
            }

            if (app.m_Toolbar["genemode"] == true) tbar.push(geneModeButton);
            if (app.m_Toolbar["slimmode"] == true) tbar.push(slimModeButton);

            tbar.push('->',
                {hidden: app.m_NeedAlignmentView == false, text: 'Sync Alignment View', scope: this,
                 handler: function(){app.onSyncAlignView(this);} },
                {iconCls: 'xsv-full_view', tooltip:'Open Entrez View', scope: this, hidden: (app.m_Embedded === true || app.m_NoViewHeader !== true),
                 handler: function(b,e){this.openFullView();}}
            );

            if (app.m_Toolbar['tools'] == true) {
                var menu_items = [
                    {text: 'Go To', iconCls: 'xsv-goto_position', tooltip: 'Go To Position/Range',ariaLabel: 'Go To Position/Range', scope: this,
                     handler: function(){ this.gotoPositionDlg(); this.pingClick('2-2-7-0');}},
                    {text: 'Search', iconCls: 'xsv-search-button', tooltip: 'Search features, components or sequences',ariaLabel: 'Search features, components or sequences', scope: this,
                     handler: function(){ app.showSearchParamsDlg(this); this.pingClick('2-2-7-1');}},
                    {text: 'Flip Strands', iconCls: 'xsv-flip-strands', tooltip:'Flip Sequence Strands',ariaLabel:'Flip Sequence Strands',
                     pressed: this.getFlip(), enableToggle:true, hidden: !this.canFlip(), scope: this,
                     handler: function(){ this.flipStrand(); this.pingClick('2-2-7-2');}},
                    {text: 'Markers', iconCls:'xsv-markers', scope: this,
                     handler: function(){ app.showMarkersDlg(this); this.pingClick('2-2-7-3');}},
                    {text: 'Set Origin', iconCls:'xsv-origin', scope:this,
                     handler: function(){ app.showOriginDlg(null); this.pingClick('2-2-7-4');}},
                    {text: 'Sequence Text View', iconCls:'xsv-new_fasta', tooltip: 'Create New Sequence Text View',ariaLabel: 'Create New Sequence Text View', scope: this,
                     handler: function(){ app.createTextView(this); this.pingClick('2-2-7-5');}}
                ];

                if (app.mf_MultiPanel) {
                    menu_items.push(
                        '-',
                        {text:'Add new Panel', iconCls:'xsv-new_view', tooltip:'Create New Graphical Panel', ariaLabel:'Create New Graphical Panel',
                            scope:this, handler:function() { this.addNewPanel('2-2-7-6'); }
                        },

                        {text:'Add new Panel on Range', iconCls:'xsv-new_view', tooltip:'Create New Graphical Panel on Selected Range',
                            ariaLabel:'Create New Graphical Panel on Selected Range', itemId: 'addRPanel',
                            scope:this, handler:function() { this.addNewPanel('2-2-7-7', this.getTotalSelectedRange()); }
                        },
                        {text: 'Show/Hide Status bar', scope: this, handler: function() {
                            var sb = this.m_View.down('#statusBar');
                            sb[sb.isVisible() ? 'hide' : 'show']();
//                            localStorage.setItem('NCBI/SV/statusBar', sb.isVisible() ? 'hide' : 'show');
//                            sb[localStorage.getItem('NCBI/SV/statusBar')]();
                        }}
                    );
                }
                var noPBlast = (seqInfo.acc_type !== 'DNA');
                menu_items.push(
                    '-',
                    {text:'BLAST and Primer Search', menu:new Ext.menu.Menu({items:[

                    {text:'BLAST Search (Visible Range)', iconCls:'xsv-blast', scope: this,
                     handler: function() {
                        app.blast(this, [this.m_VisFromSeq, this.m_VisFromSeq + this.m_VisLenSeq - 1]);
                        this.pingClick('2-2-7-8-0');}},
                    {text:'Primer BLAST (Visible Range)', iconCls:'xsv-primer', scope: this, disabled: noPBlast,
                     handler: function() {
                        this.primerBlast(false, [ [this.m_VisFromSeq, this.m_VisFromSeq + this.m_VisLenSeq - 1] ]);
                        this.pingClick('2-2-7-8-1');}},
                    {text:'BLAST Search (Selection)', iconCls:'xsv-blast', scope: this,
                     handler: function() {this.blastSelection(); this.pingClick('2-2-7-8-2');}},
                    {text:'Primer BLAST (Selection)', iconCls:'xsv-primer', scope: this, disabled: noPBlast,
                     handler: function() {this.primerBlast(); this.pingClick('2-2-7-8-3');}}
                    ]})},

                    '-',

                    {text: 'Preferences', iconCls: 'xsv-configure', scope: this,
                        handler: function() {this.pingClick('2-2-7-11'); this.setPreferences();}}
                );

                var tools_menu = new Ext.menu.Menu({ items: menu_items, listeners: { scope: this,
                    beforeshow: function(menu) { 
                        if (this.m_App.mf_MultiPanel)
                            menu.down('#addRPanel')[this.m_RangeSelectionSet.length ? 'enable' : 'disable']();
                    }
                }});

                tbar.push({
                    text: 'Tools', iconCls: 'xsv-tools', tooltip: 'Tools',ariaLabel : 'Tools',itemId: 'btnTools', scope: this,
                    handler: function() {this.pingClick('2-2-7');},
                    menu: tools_menu
                 });
            }

            var cfgButt = {
                text: 'Tracks', iconCls: 'xsv-config', tooltip:'Configure tracks',ariaLabel:'Configure tracks',
                itemId: 'btnTracks', scope: this};
            app.m_Toolbar['trackSets'] = (!seqInfo.disable_track_sets && TMS.TrackSets != undefined && app.m_AppContext != null)
            if (app.m_AppContext == 'SV-standalone' && !app.m_AssmContext) app.m_Toolbar['trackSets'] = false;
            var menu_items = [
                {text: 'Configure Tracks', tooltip: 'Configure current set of tracks',
                 ariaLabel : 'Configure current set of tracks',scope: this,
                 iconCls: 'xsv-config', handler: function(e) { app.showTracksConfigDlg('2-2-8-0'); }},
                {text: 'Configure Track Hubs', tooltip: 'Configure Remote Track Hubs and Tracks',
                 disabled: !seqInfo.top_level,
                 ariaLabel: 'Configure Remote Track Hubs and Tracks', scope: this,
                 handler: function () {app.showTrackHubDlg(); this.pingClick('2-2-8-TH');}}];

            if (app.m_NoConfDlg || app.m_Toolbar["config"] != true) menu_items.shift();
            if (app.m_Toolbar['thub'] != true)  menu_items.pop();

            if (app.m_Toolbar['trackSets'] == true) {
                var icon = 'xsv-search-loading';
                cfgButt.handler = function(butt) {
                    if (butt.menu.hidden) return;
                    this.pingClick('2-2-8');
                    if (!SeqView.requestTrackSets) {
                        delete app.m_currentTrackSetId;
                        SeqView.configureTrackSets(app);
                    }
                    var cmp = butt.down('#tsNCBI').disable();
                    cmp.setIconCls(icon);
                    cmp.menu.removeAll();
                    cmp = butt.down('#tsMy').disable();
                    cmp.setIconCls(icon);
                    cmp.menu.removeAll();
                    SeqView.requestTrackSets(this.processTrackSets.bind(this));
                }
                if (menu_items.length) menu_items.push('-');
                menu_items.push({text: 'NCBI Recommended Track Sets', itemId: 'tsNCBI',
                                iconCls: icon, menu: new Ext.menu.Menu({items: []}), disabled: true},
                                {text: 'My NCBI Track Collections', itemId: 'tsMy',
                                iconCls: icon, menu: new Ext.menu.Menu({items: []}), disabled: true},
                                {text: 'FAQ', scope: this,
                                iconCls: 'xsv-question', handler: function() {
                                    this.pingClick('2-2-8-3', 'goto_FAQ_page');
                                    window.open(SeqView.webNCBI + 'tools/sviewer/faq/#tracksets');}});
                cfgButt.menuAlign = 'tr-br?';
            }
            if (menu_items.length) {
                cfgButt.menu = new Ext.menu.Menu({ items: menu_items });
                tbar.push('-', cfgButt);
            }
            if (app.m_Toolbar['download'] == true) {
                var menu_items = [
                    {text:'Download FASTA', iconCls:'xsv-zoom_seq',
                        menu:new Ext.menu.Menu({items:[
                            {text:'FASTA (Visible Range)', scope: this,
                                handler: function () { this.pingClick('2-2-9-visFA'); this.downloadData(false, "fasta", null); }},
                            {text: 'FASTA (All Markers)', scope: this,
                                handler: function () { this.pingClick('2-2-9-markFA'); this.downloadAllMarkers("fasta"); }},
                            {text:'FASTA (All Selections)', scope: this,
                                handler:function() {this.pingClick('2-2-9-selFA'); this.downloadData(true, "fasta", null);}}
                        ]})},
                    {text:'Download GenBank Flat File', iconCls:'xsv-file', hidden: app.m_ViewParams.acc_type == 'protein',
                        menu:new Ext.menu.Menu({items:[
                            {text: 'Visible Range', scope: this,
                                handler: function() {this.pingClick('2-2-9-visFF'); this.downloadData(false, "flat", null);}},
                            {text:'All Selections', scope: this,
                                handler:function() {this.pingClick('2-2-9-selFF'); this.downloadData(true, "flat", null);}}
                        ]})},
                    {text:'Download Track Data', icon: SeqView.base_url + 'images/grid.png', hidden: app.m_noTrackData, scope: this,
                        handler: function() {
                            var tracks = this.m_App.getActiveTracks().filter(SeqView.TM.isDownloadable, this);
                            this.downloadTblFile(tracks);
                        }},
                    {text:'Printer-Friendly PDF/SVG', iconCls:'xsv-printer', disabled: app.m_NoPDF, scope: this,
                        handler: this.downloadImgFile}
                ];
                var downloads_menu = new Ext.menu.Menu({ items: menu_items, listeners: { scope: this,
                    beforeshow: function(menu) {

                    }
                }});

                tbar.push({
                    text: 'Download', iconCls: 'xsv-download', tooltip: 'Download', ariaLabel : 'Download', itemId: 'btnDownload', scope: this,
                    handler: function() {this.pingClick('2-2-12');},
                    menu: downloads_menu
                });
            }
            tbar.push(reloadButton);

            var first_graphic = true;
            for (var i = 0; i < app.m_Views.length; i++) {
                var view = app.m_Views[i];
                if (view && view != this && view.isGraphic()) {
                    first_graphic = false;
                    break;
                }
            }
            if (first_graphic) {
                if( app.m_Toolbar["help"] == true ){
                    tbar.push(
                        {iconCls: 'xsv-question',
                         tooltip:'Help',
                         ariaLabel:'Help',
                         layout: {type:'vbox'},
                         scope: this,
                         handler: function() {this.pingClick('2-2-10');},
                         menu: new Ext.menu.Menu(
                            {defaultOffsets: [-70,0], scope: this,
                             items: [
                                {text: 'Help', iconCls: 'xsv-question', scope:this,
                                 handler: function() {SeqView.showHelpDlg(); this.pingClick('2-2-10-1');}},
                                {text: 'Legends', icon_Cls: 'xsv-gene_mode', scope: this,
                                 handler: function() {SeqView.showHelpDlg('legends'); this.pingClick('2-2-10-4');}},
                                {text: 'Link to View', iconCls: 'xsv-link_to_page', scope: this,
                                 handler: function() {app.showLinkURLDlg(); this.pingClick('2-2-10-3');}},
                                {text: 'Feedback', iconCls: 'xsv-feedback', scope: app,
                                 handler: function() {this.showFeedbackDlg('2-2-10-0');}},
                                {text: 'About', scope: this,
                                 handler: function() {SeqView.showAboutMessage(this.m_Idx); this.pingClick('2-2-10-2');}}
                            ]})
                        }
                    );
                }
            } else {
                tbar.push(
                    {iconCls: 'tb-close', tooltip:'Close', ariaLabel:'Close', scope: this, handler:function() { this.remove(); }}
                );
            }
            var tools = [];
            if (app.m_Embedded === false) {
                tools.push(this.createCloseTBar());
            } else {
                tools.push({id:'help',qtip: 'Help', handler: function(){this.pingClick('2-2-10-1'); SeqView.showHelpDlg();}})
            }
            if (this.m_View) {
                this.m_View.destroyDockedItems();
                this.m_View.addDocked({xtype:'toolbar', dock:'top', items:tbar});
            } else {
                this.m_View = app.addView({
                    collapsible: true, title: 'Loading...', header: false,
                    tools: tools, view: this, tbar:tbar,
                    html: {
                        tag: 'div', id: this.m_DivId,
                        cls: 'graphical_div sv-drag sv-highlight sv-dblclick'
                    }
                })
            }
            this.m_View.addDocked({
                xtype: 'toolbar',  dock: 'bottom', itemId: 'statusBar',
                    items: [
                        {itemId: 'status', height:16, xtype: 'tbtext', text: ''},
                        {itemId: 'range', height:16,xtype: 'tbtext', text: ''},
                        '->',
                        {itemId: 'feedback', xtype: 'button', text: '', scope: app, iconCls: 'xsv-feedback', tooltip: 'Feedback',
                         handler: function() { this.showFeedbackDlg('10-2'); }},
                        {itemId: 'tracks', xtype: 'button', text: '', scope: app, iconCls: 'xsv-config',
                         handler: function (b, e) { SeqView.debugMode = e.ctrlKey; this.showTracksConfigDlg('10-1'); delete SeqView.debugMode;}}]}
            );
        }
        if (!this.bodyDiv) {
            this.bodyDiv = Ext.get(this.m_DivId);
            this.m_vLine = this.bodyDiv.appendChild(new Ext.Element(document.createElement('div')));
            this.m_posTT = new Ext.tip.ToolTip({target: this.bodyDiv, dismissDelay: 0, hidden: true, html: 'Loading...'
//            ,listeners: {show: function(){ 
//                    console.log(this.pageX, this.x);
//            }, sco_pe: this}
//                ,disabled: true
            });
            this.m_vLine.set({id: 'lineRuler' + this.m_Idx,
                 style: 'position:absolute;border:1px solid lightgray;border-right-width:0;top:0;width:0;z-index:10;display:none;'});
            this.bodyDiv.on({ scope: this,
                pointerdown:this.onPointerDown,
                pointerup:  this.onPointerUp,
                dblclick:   this.onDblClick,
                pointermove:this.onPointerMove,
                contextmenu:this.onContextMenu,
                pointerleave: function() {this.m_vLine.setStyle('display', 'none');}
            });
        }
        this.updateMolPlacement();
        this.updateStatusBar('feedback', '');
        this.startImageLoading(this.m_FromSeq, this.m_LenSeq, {fire_event:false});
//        this.updateStatusBar('link', '<a href=\"https://www.ncbi.nlm.nih.gov/tools/sviewer/\" target=\"_blank\">Help</a>');
    },

    addNewPanel: function(pingID, range) {
        var view = new SeqView.Graphic(this.m_App);
        if (range) {
            this.removeRangeSelection();
            view.m_FromSeq = range[0];
            view.m_LenSeq = range[1] - range[0] + 1;
        }
        this.m_App.registerView(view);
        if (pingID) this.pingClick(pingID);
    },

    updateMolPlacement: function() {
        var seqInfo = this.m_App.m_Config.SeqInfo || {};

        if (!seqInfo.placements || !seqInfo.placements.length) return;
        var ttl = this.m_View.down("#btnTitle");
        if (!ttl) return;

        var menu = ttl.getMenu(),
            el = menu.down('#molPlacement');
        if (el) menu.remove(el);
        else menu.add('-');

        Ext.each(seqInfo.placements, function() {
            var link = SeqView.base_url + '?id=' + this.accession + '&v=' + (this.from + 1) + ':' + (this.to + 1);
            menu.add({
                text: this.label || this.accession + ': ' + (this.from + 1) + '..' + (this.to + 1),
                tooltip: this.help,
                handler: function() { window.open(link); }
            });
        });
    },

    updateStatusBar: function(prop, text, icon) {
        var st = this.m_View.down('#' + prop);
        if (st) {
            if (typeof text == 'string') st.setText(text);
            if (typeof icon == 'string') st.setStyle('background', icon + (icon ? ' no-repeat' : ''));
            st.setStyle('padding', (st.el.getStyle('background-image') != 'none' ? '0 22px' : '0 3px'));
        }
        return st || {hide: Ext.emptyFn, show: Ext.emptyFn};
    },

    updateStatusInfo: function() {
        var warnTT = this.m_App.m_Config.checkWarnings() || '';
        if (warnTT) {
            var btn = this.m_View.down('#assmInfo');
            if (btn) btn.setTooltip('WARNING!<br>' + warnTT);
        }

        var warn = (this.m_App.m_appWarning || '') + warnTT;
        delete this.m_App.m_appWarning;

        var st = this.updateStatusBar('status', warn ? 'Warning' : '',
            warn ? ('url(' + SeqView.ExtIconLoading.slice(5, -18) + 'shared/warning.gif)') : '');
        if (warn) {
            if (!st.tooltip) st.tooltip = Ext.create('Ext.tip.ToolTip', {target: st.id, title: 'Warning'});
            st.tooltip.setHtml(warn);
        }
        else if (st.tooltip) st.tooltip = st.tooltip.destroy();
        this.updateStatusBar('tracks', 'Tracks shown: ' + this.m_graphTracks.length + '/' + this.m_App.m_Config.TrackConfig.length).show();
    },

    processTrackSets: function(sets) {
        var scope = this;
        var prepareTT = function() {
            var tt = this.tt;
            if (!tt) {
                if (this.type == 'default') {
                    this.tracks = scope.m_App.m_defaultTrackSet || [];
                }
                Ext.each(this.tracks, function(t, i) {
                    tt += '<p style="margin:0;' + (i%2 ? '">' : 'background-color:gainsboro;">') + t.GetName() + '</p>';
                });
            }
            if (tt) this.tooltip = new Ext.ToolTip({target: this.el, title: this.text, html: tt,
                anchor: 'right', margin: 0, track_Mouse: true, dismissDelay: 0});
        }
        var usrSets = [];
        var tmsSets = [{text: 'Reset to Default Tracks',
                        type: 'default', tsId: 'defaultTS', tt: '',
                        handler: function() { scope.applyTrackSet(this); },
                        listeners: { afterrender: prepareTT }}];
        Ext.each(sets, function() {
            var tSet = {text: this.GetName(),
                        tracks: this.GetTracks(), tt: '', tsId: this.GetId(),
                        handler: function() { scope.applyTrackSet(this); },
                        listeners: { afterrender: prepareTT }};
            switch (this.GetTrackSetType()) {
                case TMS.TrackSetType.TMS: tmsSets.push(tSet); tSet.type = 'tms'; break;
                case TMS.TrackSetType.myNCBI_Collection: usrSets.push(tSet); tSet.type = 'MyNCBI'; break;
            }
        });

        var cmp = this.m_View.down('#tsNCBI').enable();
        var addMenuItem = function() {
            if (scope.m_App.m_currentTrackSetId == this.tsId) this.iconCls = 'xsv-search-results';
            cmp.menu.add(this);
        }
        cmp.setIconCls('xsv-seq-logo');
        Ext.each(tmsSets, addMenuItem);
        cmp = this.m_View.down('#tsMy').enable();
        cmp.setIconCls('xsv-seq-logo');
        if (usrSets.length) {
            Ext.each(usrSets, addMenuItem);
            cmp.menu.add('-');
        }
        var inMyNCBI = MyNCBI.User.IsLoggedIn();
        cmp.menu.add({text: !inMyNCBI ? 'Sign in to NCBI' : 'Manage collections', scope: this,
                handler: function(){
                    this.pingClick('2-2-8-2', inMyNCBI ? 'NCBI_sign_in' : 'collections_manage');
                    window.open(SeqView.webNCBI + 'myncbi/collections');}});
        if (inMyNCBI)
            cmp.menu.add({text: 'Save current tracks', disabled: this.m_App.m_noSaveTracks, handler: function() { scope.saveTrackSet(); }});
    },

    saveTrackSet: function() {
        var idList = [],
            dOpts = {};

        this.pingClick('2-2-8-1', 'saveTrackSet');
        Ext.each(this.m_App.getActiveTracks(), function() {
            dOpts[this.id] = SeqView.TM.getTrackDisplayOptions(this).slice(0);
            idList.push(this.id);
        });
        if (!idList.length) return;
        TMS.GetTracksById(idList).done(function(){
            var tracklist = this;
            Ext.MessageBox.prompt('Save Track Collection', 'Please enter collection name:', function(btn, name) {
                if (btn!='ok'  || name.length==0) return;
                var tracks = [];
                Ext.each(tracklist.GetTracks(), function() {
                    var dti = new TMS.DisplayTrackInfo(this);
                    dti.SetDisplayOptions(dOpts[this.GetTMSId()]);
                    tracks.push(dti);
                });

                TMS.TrackSets.TrackSetService.CreateTrackset(name, tracks).done(function(trackset_id) {
                    console.log('trackset ' + trackset_id + 'has been added');
                });

            });
       });
    },

    applyTrackSet: function(tSet) {
        this.ping({'jsevent': 'click', 'sv-area': '2-2-8-' + tSet.type,
                   'sv-event': 'applyTrackSet', trackSetID: tSet.tsId});
        var trx = '',
            ts = {tracks: [], ignoreURL: true};
        Ext.each(tSet.tracks, function() {
            ts.tracks.push({displayOptions: this.GetDisplayOptions(), id: this.GetTMSId()});
        });
        this.m_App.fireEvent('trackSet_requested', ts);
        Ext.each(ts.tracks, function() {
            trx += '[id:' + this.id + (this.displayOptions ? ',' + this.displayOptions : '') + ']';
        });
        this.loadTrackSet(trx || '[amend]', tSet.tsId, ts.ignoreURL);
    },
    
    updateTracks: function(tracks) {
        if (this.deferredUpdate) clearTimeout(this.deferredUpdate);
        if (this.m_Loading) {
            this.deferredUpdate = Ext.defer(this.updateTracks, 100, this, [tracks]);
            return;
        }
        delete this.deferredUpdate;
        this.loadTrackSet(tracks);
    },

    loadTrackSet: function(tracks, tsId, ignoreURL) {
        var app = this.m_App;
        this.setGeneMode(false);
        var dflt = app.defaultConfig.m_TracksFromURL || '';
        if (ignoreURL) dflt = '';

        if (tsId !== 'defaultTS')
            app.m_TracksFromURL = dflt.replace(/]/g, ',shown:false]') + tracks;
        else
            app.m_TracksFromURL = tracks + dflt;
        app.m_currentTrackSetId = tsId;
        this.m_View.down('#btnReload').disable();

        app.m_Config.load({forcereload: true,
            callback: { scope: this,
                success: function(res){ app.fireEvent('configuration_changed', app); this.refresh(); },
                failure: function(data, text, res) { app.infoMissed(data, text, res); }
            }
        });
    },

    clear: function() {
        this.m_History = [];
        this.m_HistIx = 0;
        this.mf_HistoricalBack = false;

        this.m_Selection = null;
        this.m_RangeSelectionSet = [];
        this.m_SelectedSig = null;

        if( this.m_App){
            this.m_App.fireEvent( 'selection_changed', this );
        }

        this.m_ScrollPixBeforeLoad = 0;

        SeqView.Graphic.superclass.clear.call(this);
    },

    refreshCTRL: function(options) {
        if (false) {
            this.m_App.GI = 'NC_' + ('00000' + (Math.random()*22 + 1).toFixed()).slice(-6);
            delete this.m_App.m_Config.TrackConfig;
            
            this.m_App.m_Config.load({forcereload: true, notrackconfig: 1,
                callback: { scope: this,
                    success: function() {
                        var seqInfo = this.m_App.m_Config.SeqInfo;
                        var trackInfo = this.m_App.getActiveTracks();
                        var appParams = this.m_App.getAppParams();
                        this.m_App.loadSequence(seqInfo, trackInfo, appParams);
                    },
                    failure: function() {this.reload();}
                }
            });
            return;
        }
        //compare providers and affinity_tokens
        var count = 0;
        Ext.each(this.m_App.m_Config.TrackConfig, function(t) {
            t.dp = t.data_provider;
            t.at = t.affinity_token;
            delete t.data_provider;
            delete t.affinity_token;
        });
        Ext.each(this.m_App.m_Config.TrackConfig, function(t) {
            if (t.dp != this.setProviderAndAffinity(t))
                console.log(count++, t.dp, t);
            if (t.at != t.affinity_token)
                console.log(count++, t.at, t);
            return count < 10;
        }, this.m_App.m_Config);
        if (!count) console.log('No errors detected');
    },

    refresh: function(options) {
        if (this.m_VisLenSeq == 0)
            return;
        if (options && options.buttonAction) {
            delete options;
            if (!this.m_Loading && this.m_graphTracks) {
                Ext.each(this.m_graphTracks, function(t) { delete t.tParams;});
            }
        }
        this.startImageLoading(this.m_VisFromSeq, this.m_VisLenSeq, options);
        if (options && options.resize) this.updateTitle();
    },

    // remove hover selection decoration
    removeFloatingSelection: function() {
       if (this.m_Selection) {
           this.m_Selection.remove();
           this.m_Selection = null;
       }
    },

    selectRange: function(sel_range) { this.selectRangeSet([sel_range]); },

    selectRangeSet: function(range_set, noEvent) {
        if (!range_set || !range_set.length) return;
        this.removeRangeSelection(null, noEvent);
        for (var x in range_set) {
            this.m_RangeSelectionSet.push(new SeqView.RangeSelection(this, range_set[x]));
        }
    },
    //removes RangeSelections and SelectedRangeSets when no pinned range tooltip associated with them
    removeSelectionsWithNoPinnedToolTips: function() {
        var length = this.m_RangeSelectionSet.length;
        for (var x = length - 1; x >= 0; x--) {
            var RSSx = this.m_RangeSelectionSet[x];
            RSSx.remove();
            this.m_RangeSelectionSet.splice(x, 1);
        }
        if (length != this.m_RangeSelectionSet.length) this.m_App.fireEvent('selection_changed', this);
    },

    popRangeElement: function(range_elt) {
        var sel_set = this.m_RangeSelectionSet;
        for (var x in sel_set) {
            if (sel_set[x] == range_elt) {
                sel_set.splice(x, 1);
                break;
            }
        }
    },

    removeRangeSelection: function(range, noEvent) {
        var sel_set = this.m_RangeSelectionSet;
        if (range) {
            for (var x in sel_set) {
                if (sel_set[x].range != range) continue;
                sel_set.splice(x, 1)[0].remove();
                break;
            }
        } else {
            Ext.each(sel_set, function(rs) { rs.remove(); });
            this.m_RangeSelectionSet = [];
        }
        if (!noEvent) this.m_App.fireEvent('selection_changed', this);
    },

    removeFeatMarks: function() {
        if (this.m_HighLightedFeature) {
            Ext.each(this.m_HighLightedFeature,function(sf) {
                sf.remove();
            });
        }
        delete this.m_HighLightedFeature;
    },

    addRangeSelection: function(sel_element, range, pixel_width) {
        this.m_InDragAction = false;
        if (pixel_width < 5) {
            sel_element.remove();
            return false;
        }
        sel_element.updateCoords();
        this.m_RangeSelectionSet.push(sel_element);
        this.m_App.fireEvent( 'user_changed_selection', this );
        return true;
    },

    getTotalSelectedRange: function() {
        if (!this.m_RangeSelectionSet.length) return null;
        var min_pos = this.m_App.m_SeqLength;
        var max_pos = 0;
        Ext.each(this.m_RangeSelectionSet, function() {
            var rng = this.range;
            max_pos = Math.max(max_pos, rng[0], rng[1]);
            min_pos = Math.min(min_pos, rng[0], rng[1]);
        });
        return [min_pos, max_pos];
    },

    getSelectedRangeSet: function() {
        var srs = [];
        Ext.each(this.m_RangeSelectionSet, function(){ srs.push(this.range); });
        return srs;
    },

    startImageLoading: function(from, len, options) {
        if (!arguments.length) {
            if (!this.m_deferLoadParams) return;
            from = this.m_deferLoadParams[0];
            len = this.m_deferLoadParams[1];
            options = this.m_deferLoadParams[2];
        }
        if (this.m_Loading) {
            var info = { target: this, method: this.startImageLoading, args: arguments };
            this.m_App.fireEvent('defer_image_load', info);
            if(!info.ignore) {
                if (!this.m_deferLoadParams) 
                   Ext.defer(this.startImageLoading, 1000, this);
                if (arguments.length) this.m_deferLoadParams = arguments;
            }
            return;
        };
        delete this.m_deferLoadParams;
        this.clearCache();
        this.updateGraphicData(from, len, options);
    },

    // Check that the range fits sequence boundaries,
    // if not, adjust it appropriately
    adjustRequest: function(from, len) {
        if (from < 0) from = 0;
        if (from + len > this.m_App.m_SeqLength) {
            var over = from + len - this.m_App.m_SeqLength;
            from -= over;
            if (from < 0) {
                from = 0;
                len = this.m_App.m_SeqLength;
            }
        }
        return [from, len];
    },


    processGraphicParams: function(params, reqNum) {
        this.m_App.m_timeStamps.graphicDataUpdateTS = new Date().getTime();
        var rqCfg = this.rqParams ? this.rqParams.rqCfg : {};
        var reuse = true === Ext.each(['assm_context','color','content','decor','flip','id','label','len','select','spacing','width', 'hidden_settings'],
            function(prop) { return params.rqCfg[prop] == rqCfg[prop]; });
        reuse = reuse && this.m_VisFromSeq >= rqCfg.from && this.m_VisFromSeq + this.m_VisLenSeq <= rqCfg.from + rqCfg.len;
        reuse = reuse && !this.noParallel;
        this.m_showGeneModeButton = false;
        var newTrx = [];
        this.m_graphTracks = this.m_graphTracks || [];
        this.m_FromCgi = {areas:[], images:[]};
        this.updateStatusBar('status', 'Loading...', SeqView.ExtIconLoading);
        this.updateStatusBar('tracks').hide();

        Ext.each(this.m_App.getActiveTracks(), function(nTrk, idx) {
            if (!nTrk.shown) return;
            this.m_showGeneModeButton |= nTrk.category.name == 'Genes';
            nTrk.hasOpts = (nTrk.choice_list || nTrk.check_boxes || nTrk.text_boxes || nTrk.range_control) != undefined;
            if (nTrk.labels) {
                Ext.each(nTrk.labels, function() { Ext.removeNode(this); });
                delete nTrk.labels;
            }

            this.m_graphTracks.forEach(function(t, i) {
                if (t.id == nTrk.id && t.signature == nTrk.name) {
                    delete t.status;
                    t.cfgIdx = nTrk.cfgIdx;

                    nTrk = t;
                    delete this.m_graphTracks[i];
                    return false;
                }
            }, this);
            if (!reuse || !nTrk.tParams || nTrk.tParams != params.tracks[idx]) {
                Ext.each(nTrk.labels, function() { Ext.removeNode(this); });
                delete nTrk.labels; 
                delete nTrk.image;
            }
            nTrk.tParams = params.tracks[idx];
            newTrx.push(nTrk);
        },this);
        this.removeGraphTracks();
        this.m_graphTracks = newTrx;

        if (reuse) { 
            this.rqParams.rqCfg.origin = params.rqCfg.origin || 0;
            params.rqCfg = this.rqParams.rqCfg;
        } else delete this.m_seqPrimers;
        this.rqParams = params;
        this.tmpScrWidth = this.getBodyWidth();
        this.m_LenSeq = params.rqCfg.len;
        this.m_Width = params.rqCfg.width;
        this.m_bpPix = (this.m_Width || SeqView.ChunkWidth)/this.m_LenSeq;
        this.m_FromSeq = params.rqCfg.from;
        this.m_Height = this.m_bottomRulerHeight = this.m_topRulerHeight = this.m_topRulerHeight || 18;
        this.clearFeatMarkers();
        this.updateTrackData({tParams : '[key:ruler]'}, -1); //request ruler
        delete rqCfg.origin;
        Ext.each(newTrx, this.updateTrackData, this);
        this.updateViewHeight(-1);
    },


    updateScrollPos: function() {
        this.m_ScrollPix -= this.m_ScrollPixBeforeLoad ? this.m_ScrollPixBeforeLoad : this.m_ScrollPix; 
        this.m_ScrollPixBeforeLoad = 0;
        // calculating initial offset
        this.m_ScrollPix -= this.getFlip()
            ? Math.max(0, this.seq2Pix(this.m_ReqFrom - 0.5) - this.getBodyWidth())
            : Math.min(Math.max(0, this.m_Width - this.getBodyWidth()), this.seq2Pix(this.m_ReqFrom - 0.5));
    },

    adjustViewWidth: function() {
        var scrWidth = this.getBodyWidth();
        if (scrWidth != this.tmpScrWidth)
            this.m_graphTracks.forEach(function(t) {
                if (t.tnEl) {
                    if (t.tnEl.getWidth() == scrWidth) return;
                    t.tnEl.setWidth(scrWidth);
                    var btns = t.tnEl.query('button[class=xsv-jumpbutton]', false);
                    if (btns.length == 2) {
                        var divWidth = (scrWidth >> 1) - 26; 
                        var div = t.tnEl.query('div',false)[0];
                        if (div) div.setWidth(divWidth);
                        btns[0].setLeft(divWidth);
                        btns[1].setLeft(divWidth + 32);
                    }
                }
                if (t.labels) t.labels[1].setRight(12 - scrWidth);  
            });
        delete this.tmpScrWidth;
    },

    updateViewHeight: function(startTrk) {
        var gTrx = this.m_graphTracks,
            height = this.m_topRulerHeight + this.m_bottomRulerHeight;
        if (startTrk !== null) {
            var track = gTrx[startTrk];
            lastBottom = track ? track.div.getBottom(true) : this.m_topRulerHeight;
            
            for (var i = startTrk + 1; i < gTrx.length; i++) {
                var trk = gTrx[i];
                if (!trk || !trk.div || !trk.div.dom) continue;
                trk.divT = lastBottom;
                this.shiftGraphTrack(trk);
                lastBottom = trk.div.getBottom(true);
            }
        }

        if (gTrx.length) height = gTrx[gTrx.length - 1].div.getBottom(true) + this.m_bottomRulerHeight;
        this.m_Height = this.bodyDiv.getHeight();

        if (this.m_Height < height || !this.m_Loading) {
            this.bodyDiv.setHeight(this.m_Height = height);
            this.m_View.updateLayout();
        }
        if (this.bottomRuler) this.bottomRuler.setTop(this.m_Height - this.m_bottomRulerHeight);
    },

    updateTrackData: function(trk, idx) {
        if (trk.image) {
            this.tmpScrWidth = -1;
            return;
        }
        var prev,
            gTrx = this.m_graphTracks,
            params = this.rqParams,
            rqCfg = params.rqCfg,
            affTokens = params.affTokens || [];
            

        if (idx >= 0) { //not ruler so adding a placeholder
            for (var i = idx; --i >= 0;) if (prev = gTrx[i]) break; //find any of previous
            var top = prev ? prev.div.getBottom(true) : this.m_topRulerHeight;
            var div = trk.div;
            if (!div) {
                div = new Ext.Element(document.createElement('div'));
                div.applyStyles('top: ' + top + 'px; width: ' + rqCfg.view_width
                    + 'px; height: 18px; position:absolute; font-family:Monospace,Courier New,serif; font-size: 12px;'
                    + 'overflow:hidden;white-space: nowrap; text-overflow: ellipsis;'
                    + 'background:white;background:linear-gradient(#eeeeee, white, white);');
                div.setHtml(trk.display_name + ', Loading...');

                if (prev) div.insertAfter(prev.div);
                else this.bodyDiv.insertFirst(div, true);

                gTrx[idx] = { div: div, tnEl: div, display_name: trk.display_name, divT: top,
                    imgH: div.getHeight(), id: trk.id, cfgIdx: trk.cfgIdx, tParams: trk.tParams };
                div.dom.cfgIdx = trk.cfgIdx;
            } else {
                this.shiftGraphTrack(trk, top - trk.divT);
                if (trk.tnEl) {
                    var tnDiv = trk.tnEl.first() || trk.tnEl;
                    if (tnDiv.dom.tagName != 'DIV') tnDiv = trk.tnEl;
                    tnDiv.setHtml(tnDiv.getHtml().replace(trk.display_name, trk.display_name + ', Loading...'));
                }
            }
            if (!trk.tnEl || !trk.tnEl.query('button[class=xsv-hidetrack]').length) this.addXButton(div.dom);
            for (var i = idx + 1; i < gTrx.length; i++) {
                var t = gTrx[i];
                if (!t || !t.div) break;
                this.shiftGraphTrack(t, gTrx[i - 1].div.getBottom(true) - t.divT);
            }
        }

        var maxRQs = 1000;
        var submitRequest = function(trk, idx) {
            if (!trk) {
                idx = this.deferredRQs.shift();
                if (!this.deferredRQs.length) delete this.deferredRQs;
                trk = gTrx[idx];
            }
            if (!trk.tParams) return;
            rqCfg.tracks = trk.tParams;
            if (affTokens[idx]) rqCfg.affinity_token = affTokens[idx]
            else delete rqCfg.affinity_token;;
            if (trk.selectedSig && !this.noParallel) rqCfg.select = trk.selectedSig;
//rqCfg.delaytest=10;
            var rqID = this.m_App.AjaxRequestExt({
                url: params.tURLs[idx] || this.m_App.m_CGIs.Graphic, params: rqCfg, scope: this,
                callback: function(opts, success, res) { processResponse(res, idx, opts, success); }
            }, this, idx).id;
            if (rqID >= 0) this.m_Loading++;
//console.log('Submitted: ', this.m_Loading, 'Deferred: ', this.deferredRQs);
            if (idx >= 0) gTrx[idx].rqID = rqID;
            delete rqCfg.tracks;
            if (!this.noParallel) delete rqCfg.select;
        }.bind(this);

        var processResponse = function(res, idx, opts, success) {
            if (!success) this.m_App.watchdogReport(res);
            var data = success ? Ext.decode(res.responseText) : res;
            if (idx >= 0 && (!gTrx[idx] || gTrx[idx].rqID != res.requestId || !gTrx[idx].imgH)) return;

            if (!data.job_status) {
if (this.m_App.m_hiddenOptions.indexOf('ignore_areas') == -1)
                if (!data.areas || !data.areas.length) {
                    var msg = 'Incomplete data returned by back-end service',
                        app = this.m_App;
                    this.m_App.submitError('Missing or empty "areas" property returned', opts.url, res.request.requestOptions.data, '');
                    if (this.noParallel) {
                        Ext.Msg.show({ title: 'Error', message: msg, buttons: Ext.Msg.YESNO, buttonText: { yes: 'OK', no: 'Feedback'}, icon: Ext.Msg.ERROR,
                            fn: function(btn) { if (btn === 'no') app.showFeedbackDlg('2-Feedback', msg); }
                        });
                    }
                    this.addFailedTrack(idx, msg);
                    return;
                }
//                data.areas.forEach(function(a) { if (!(a.type & SeqView.AreaFlags.Track)) { a.bounds.t -= 0; a.bounds.b -=0; }});
                var it  = {img: new Image(), idx: idx, h: data.img_height, data_height: data.data_height, areas: data.areas};
                if (idx >= 0) {
                    this.m_FromCgi.images[idx] = it;
                    delete gTrx[idx].timeoutID;
                    gTrx[idx].rqID = -1;
                }
                if (data.hairlines && !this.noParallel) {
                    it.hairlines = data.hairlines;
                    data.hairlines.forEach(function(hl){hl.pos = Math.floor(hl.pos);});
                }
                var img_url = data.img_url;
                if (img_url && img_url.charAt(0) == '?') {
                    img_url = this.m_App.m_CGIs.NetCache + img_url;
                }

                it.img_url = it.img.src = img_url;
                var onloadImg = function() {
                    if (this.updatingPosition) {
                        Ext.defer(onloadImg, 10, this);
                        return;
                    }
                    this.updateScrollPos();
                    if (idx != -1 && this.noParallel) {
                        var areas = it.areas,
                            trackConfig = this.m_App.m_Config.TrackConfig;
                        delete it.areas;
                        var images = [];
                        Ext.each(areas, function(a, idxA) {
                            if (!(a.type & SeqView.AreaFlags.Track)) return a.type & SeqView.AreaFlags.Sequence;
                            var shortSig = a.signature.split(';')[0];
                            var image = {areas: [a]};
                            Ext.each(gTrx, function(t, idxGT) {
                                if (trackConfig[t.cfgIdx].name != shortSig) return;
                                t.signature = shortSig;
                                image.idx = idxGT;
                                var pId = trackConfig[t.cfgIdx].id;
                                for (var i = idxA + 1; i < areas.length; i++) {
                                    var aa = areas[i];
                                    if (aa.parent_id && aa.parent_id == pId) image.areas.push(aa);
                                }
                                return false;
                            });
                            images.push(image);
                        });
                        this.m_Loading += images.length - 1;

                        Ext.each(images, function(image) {//backwards to define track heights 
                            it.areas = image.areas;
                            it.idx = image.idx;
                            var lastTop = image.areas[0].bounds.t
                            this.addTrackImage(Ext.apply(image, it));
                            it.h = lastTop; 
                        }, this, true);
                        console.log('done');
                    } else this.addTrackImage(it);
                    while (this.deferredRQs && this.m_Loading < maxRQs) {
                        submitRequest();
                    }
                    if (this.m_Loading) return;
                    if (SeqView.scrolled) {
                        Ext.defer(onloadImg, 10, this);
                        return;
                    }
                    this.applyResultsOnLoad();
                }.bind(this);
                it.img.onload = onloadImg;
                it.img.empty = (typeof it.areas == 'undefined');
            } else {
                var jst = data.job_status;
                if (jst == 'submitted' || jst == 'running' || jst == 'pending') {
                    opts.params = { job_key: data.job_id };
                    if (idx >= 0)
                        gTrx[idx].timeoutID = Ext.defer(function() {
                            if (!gTrx[idx] || !this.m_Loading) return;
                            if (!gTrx[idx].imgH) opts.params.Cancel = true;
                            gTrx[idx].rqID = opts.params.Cancel ? -1 : SeqView.App.simpleAjaxRequestExt(opts).id;
                        }, 2000, this);
                    else
                        Ext.defer(SeqView.App.simpleAjaxRequestExt, 1000, this, [opts]);
                } else {
                    this.addFailedTrack(idx, data.error_message || (res.status + ' ' + res.statusText), data.job_status);
                }
            }
        }.bind(this);

        if (this.m_Loading < maxRQs) {
            submitRequest(trk, idx);
        } else {
            this.deferredRQs = this.deferredRQs || [];
            this.deferredRQs.push(idx);
        }
    },


    updateGraphicData: function(from, len, options){
        from = Math.floor(Math.max(0, from));
        len = Math.floor(Math.min(len,(this.m_App.m_SeqLength - from)));

        var scrWidth = this.getBodyWidth(); // our panel width in pixels
        if (scrWidth <= 0 || !this.m_View.isVisible()) {
            // the view is hidden, just save the visible range.
            // the client will have to call refresh method to load the image
            this.m_VisFromSeq = from;
            this.m_VisLenSeq = len;
            this.m_Loading = 0;
            this.m_App.m_InitialLoading = false;
            return;
        }

        var req_from = from;
        var req_len = len;

        // Basepairs can not be wider than 1/SeqView.MinBpp (24 pixels per basepair)
        var new_req_len = Math.floor(scrWidth * SeqView.MinBpp);
        if (req_len < new_req_len) {
            new_req_len = Math.min(this.m_App.m_SeqLength, new_req_len);
            options = options || {}
            options.update_locator = true;
            // we check for range validity and adjust it later
            req_from = Math.floor(req_from - (new_req_len - req_len)/2);
            req_len = new_req_len;
        }

        // Check that view is detailed, so we can't get away with fractional
        // pixels per basepair. If there are more than 4 pixels per base, bases
        // are becoming individually visible, so fractional ppb will lead to
        // uneven spacing between base letters
        var ppb = scrWidth/req_len;
        if (ppb > 4) {
            ppb = Math.floor(ppb);
            var new_req_len = Math.floor(scrWidth/ppb);
            if (new_req_len > this.m_App.m_SeqLength) {
                new_req_len = this.m_App.m_SeqLength;
                ppb = scrWidth/new_req_len;
            }
            req_len = new_req_len;
        }
        var from_len = this.adjustRequest(req_from, req_len);
        req_from = from_len[0]; req_len = from_len[1];

        var chunk_width;
        if( req_len > 2000000 ){
            chunk_width = 3*scrWidth;
        } else if( req_len > 1000000 ){
            chunk_width = 5*scrWidth;
        } else {
            chunk_width = 8*scrWidth;
        }
        chunk_width = Math.min(SeqView.ChunkWidth, chunk_width);
        // Actually chunk can be little bigger than maximum.

        var fit_len = Math.floor( chunk_width / ppb );

        from = Math.floor( req_from - (fit_len - req_len)/2 );
        len = fit_len;

        var from_len = this.adjustRequest(from, len);
        from = from_len[0]; len = from_len[1];
        chunk_width = len * ppb;

        if (!isNaN(chunk_width)) {
            chunk_width = Math.floor(chunk_width);
        }
         // added some sanity check before calling seqgraphic.cgi
        if (isNaN(from) || isNaN(len) || isNaN(chunk_width) || chunk_width > 16384) {
            // TODO we need to understand why we got in here.
            // Nothing left than just cancel loading and return.
            this.m_Loading = 0;
            return;
        }

        len = Math.max(len, Math.min(this.m_App.m_SeqLength, 100));
        if (chunk_width == scrWidth) req_len = len;
        this.m_ReqFrom = req_from;
        this.m_ReqLen  = req_len;
        
        if (options && options.doMove && len == this.m_LenSeq && chunk_width == this.m_Width) {
            delete options.doMove;
            this.moveTo(req_from, req_len, options);
        }
        else {
            if (this.m_VisFromSeq != req_from || this.m_VisLenSeq != req_len) {
                this.updateTitle(options);
            }
            this.m_View.down('#btnReload').disable();
            document.getElementById(this.m_App.m_DivId).style.cursor = 'progress';
            // remove selection decorations
            this.removeFloatingSelection();
            this.m_SelectedRangeSet = this.m_SelectedRangeSet || this.getSelectedRangeSet();
            this.removeRangeSelection();
            this.removeFeatMarks();
            var params = this.getGraphicParams(from, len, chunk_width);
    
            if (options && options.data_changed === true) params.rqCfg.data_changed = 'true';
    
            this.m_App.watchdogStart('', '', params/*, function() { this.m_View.down('#btnReload').enable(); }.bind(this)*/);
            params.options = options || {};
            this.processGraphicParams(params, ++this.m_ReqNum);
        }
        if (this.m_TbSlider) {
            var slider_range = 100;
            var vis_len = req_len || this.toSeq()[2];
            var min_bases_to_show = this.getMostRightPix() * SeqView.MinBpp;

            var ee = Math.log(this.m_App.m_SeqLength / min_bases_to_show);
            var slider_val = slider_range * Math.log(vis_len / min_bases_to_show) / ee;

            this.m_TbSlider.setValue(100 - slider_val);
        }
    },


    getGraphicParams: function(from, len, width, singleImg) {
        var app = this.m_App;
        var params = {id: app.GI, client:'seqviewer', width: width,
                      view_width: this.getBodyWidth(), from:from, len:len};
        Ext.apply(params, app.m_GraphicExtraParams);


        //if (this.m_SelectedSig) params.select = this.m_SelectedSig;
        if (this.getFlip()) params.flip = 'true';

        app.addURLParams(params);
        var actTrx = app.getActiveTracks(),
            sm = this.m_slimMode ? ',show_title:false' : '',
            trkParams = {tURLs: [app.m_CGIs.Graphic]};
        this.noParallel = (this.m_SelectedSig && this.noParallel) || app.m_hiddenOptions.indexOf('no_parallel') >= 0;
        if (!this.noParallel && this.m_SelectedSig) {
            var num = 0;
            if (this.m_graphTracks)
                this.noParallel = this.m_graphTracks.some(function(t){num += t.hasOwnProperty('selectedSig'); return num > 1;});
            else {
//                num = this.m_SelectedSig.indexOf(';');
                this.noParallel = true; //(num > 0 && num + 1 < this.m_SelectedSig.length);
            }
        }

        if (singleImg || this.noParallel) {
            if (this.m_SelectedSig) params.select = this.m_SelectedSig;
            trkParams.tracks = [SeqView.TM.tracksArrayToString(actTrx, true, singleImg, {extraStr: sm})];
            Ext.each(actTrx, function() {
                if (this.pURL == app.m_CGIs.Graphic) return;
                trkParams.tURLs = [this.pURL];
                return false;
            });
        }
        else trkParams = SeqView.TM.tracksToArrayOfStrings(actTrx, sm);
        trkParams.tIDs = [];
        Ext.each(actTrx, function(t) {
            if (this.seq_data_key) params.key = (params.key || '') + ',' + this.seq_data_key;
            trkParams.tIDs.push(t.id);
        });
        if (this.noParallel) delete trkParams.tIDs;
        if (params.key) params.key = params.key.replace(/^,/, '');         
        trkParams.tracks = trkParams.tracks || ['[key:ruler]'];

        if (app.m_ViewLabels) params.markers = app.m_ViewLabels;
        if (app.m_SnpFilter) params.snp_filter = app.m_SnpFilter;

        params.markers = app.m_MarkersInfo.getMarkersData(true, singleImg);
        return Ext.apply({rqCfg: params}, trkParams);
    },

    clearFeatMarkers: function(keep) {
        var fm = this.m_featMarkers;
        if (fm && fm.div) this.bodyDiv.removeChild(fm.div);
        if (!keep) delete this.m_featMarkers;
    },

    drawFeatMarkers: function() {
        if (this.m_graphTracks.length <= 1) return;
        if (!this.m_featMarkers) return; //this.m_SelectedSig = null;
        this.clearFeatMarkers(true);

        var div = this.m_featMarkers.div = document.createElement('div');
        var trx = this.m_graphTracks,
            top = trx[0].divT,
            height = this.m_Height - top - this.m_bottomRulerHeight,
            left = this.m_featMarkers.left;
        div.setAttribute('class', 'sv-drag sv-dblclick');
        div.setAttribute('style', 'position:absolute; top:' + top
            + 'px; left:' + (left + this.m_ScrollPix)
            + 'px; width:' + this.m_featMarkers.width
            + 'px; height:' + height + 'px;');
        var tmpl = '<div class="sv-drag sv-dblclick"'
            + ' style="border-width:1px; border-left-style:solid; position:absolute;';

        this.m_featMarkers.hls.forEach(function(hls){
            var arrHLs = hls.hairlines;
            var tHeight = hls.track.divT - top,
                bTop = hls.track.divT + hls.track.imgH - top,
                bHeight = height - bTop;

            for(var i = 0; i < arrHLs.length; i += 2) {
                 var lHL = arrHLs[i],
                     rHL = arrHLs[i + 1] || lHL,
                     lPos = lHL.pos - left,
                     width = rHL.pos - lHL.pos + 1;

                if (i == arrHLs.length - 1) width = -1;

                if (typeof lHL.top_color == 'object') lHL.top_color = (lHL.top_color.r == 255) ? 'red' : 'grey';//rgbToHex(lHL.top_color);
                if (typeof rHL.top_color == 'object') rHL.top_color = (rHL.top_color.r == 255) ? 'red' : 'grey';//rgbToHex(rHL.top_color);
                if (tHeight)
                    div.innerHTML += tmpl + 'width:' + width + 'px; border-right-style:' + (width >= 0 ? 'solid' : 'none')
                        + ';border-left-color:' + lHL.top_color + '; border-right-color:' + rHL.top_color
                        + ';opacity: 0.4;top: 0px; height:' + tHeight + 'px; left:' + lPos + 'px;"></div>';
                if (bHeight > 0)
                    div.innerHTML += tmpl + 'width:' + width + 'px; border-right-style:' + (width >= 0 ? 'solid' : 'none')
                    + '; border-left-color: darkgreen; border-right-color: darkgreen; opacity: 0.7;'
                    + 'top:' + bTop + 'px; height:' + bHeight + 'px; left:' + lPos + 'px;"></div>';
            }
        });

        var view = this;
        div.onmousemove = function(e) { view.highlightElement(new Ext.event.Event(e)); };
        this.bodyDiv.appendChild(div);
    },

    cleanGraphTrack: function(t) {
        if (!t) return;
        t.lgEl = Ext.removeNode(t.lgEl);
        Ext.each(t.labels, function() { Ext.removeNode(this); });
        delete t.labels;
        if ((t.div != t.tnEl)&&(t.div.dom)) t.div = Ext.removeNode(t.div);
        t.tnEl = Ext.removeNode(t.tnEl);
    },

    removeGraphTracks: function(idx) {
        var trx = this.m_graphTracks;
        if (!trx) return null;
        if (idx !== undefined) trx = trx.splice(idx, 1);
        else delete this.m_graphTracks;
        Ext.each(trx, this.cleanGraphTrack);
        return trx[0];
    },

    addSeqPrimers: function(track, primers) {
        var seq_t = primers[0] - 4 + track.divT - track.imgT;
        var seq_b = primers[1] - 8 + track.divT - track.imgT;
        //left2right sequence 5' label
        var prm = new Ext.Element(document.createElement('div'));
        var img = '<img data-qtitle="5 \'" data-qwidth="80" data-qtip="5 Prime end" aria-label="5 Prime end" src="' + SeqView.base_url + 'images/5prime_';
        prm.insertHtml('afterBegin', img + 'r.gif">');
        prm.applyStyles('position:absolute;z-index:10;left:0px;opacity:0.6;top:' + ((!this.m_Flip) ? seq_t : seq_b) +'px;');
        track.labels = [this.bodyDiv.appendChild(prm.dom)];
        //right2left sequence 5' label
        prm = new Ext.Element(document.createElement('div'));
        prm.insertHtml('afterBegin', img + 'l.gif">');
        prm.applyStyles('position:relative;z-index:10;right:' + (12 - this.getBodyWidth()) + 'px;opacity:0.6;top:'
                + ((!this.m_Flip) ? seq_b : seq_t) +'px;');
        track.labels.push(this.bodyDiv.appendChild(prm.dom));
    },

    addXButton: function(node) {
        var tpl = new Ext.Template('<button id="{id}" class="xsv-hidetrack xsv-btntrack" title="Hide track"></button>');
        var bttn = tpl.append(node, {id: node.id + 'xb'}, true);
        var handler = function(e) {
            e.stopPropagation();
            e.preventDefault();
            this.pingClick('2-0-1');
            var cfgIdx = e.target.parentNode.cfgIdx;
            if (cfgIdx == undefined) return;
            this.m_App.m_Config.TrackConfig[cfgIdx].shown = false;
            if (!this.m_Loading) {
                if (this.m_SelectedSig && !this.m_featMarkers)
                    SeqView.TM.Common.updateSeqViewApp(
                        SeqView.TM.processTracksInfo(this.m_App.m_Config.TrackConfig), this.m_App);
                else 
                    this.hideTrack(cfgIdx);
            } else {
                Ext.each(this.m_graphTracks, function(trk, idx) {
                    if (trk.cfgIdx != cfgIdx) return;
                    trk.div.setHeight(0);
                    trk.tnEl.setHeight(0);
                    trk.imgH = 0;
                    this.updateViewHeight(idx);
                    if (trk.rqID >= 0) {
                        if (trk.timeoutID) clearTimeout(trk.timeoutID);
                        delete trk.timeoutID;
                        trk.rqID = -1;
                        if (--this.m_Loading == 0) this.applyResultsOnLoad();
                    }
                    return false;
                }, this);
            }
        }
        bttn = bttn || Ext.get(node.id + 'xb');
        bttn.setListeners({touchstart: handler, click: handler, scope: this});
        return 15;
    },

    addCfgButton: function(node) {
        var tpl = new Ext.Template('<button id="{id}" class="xsv-cfgtrack xsv-btntrack" title="Track display options"></button>');
        var bttn = tpl.append(node, {id: node.id + 'cfgb'}, true);
        var handler = function(e) {
            e.stopPropagation();
            e.preventDefault();
            this.pingClick('2-0-Cfg');
            var cfgIdx = e.target.parentNode.cfgIdx;
            if (cfgIdx == undefined) return;
            SeqView.TM.modifyTrackDetails(this, cfgIdx);
        }
        bttn = bttn || Ext.get(node.id + 'cfgb');
        bttn.setListeners({touchstart: handler, click: handler, scope: this});
        return 15;
    },

    addDownloadButton: function(node, pos) {
        var tpl = new Ext.Template('<button id="{id}" class="xsv-dwnldtrack xsv-btntrack" title="Download Track Data"></button>');
        var bttn = tpl.append(node, {id: node.id + 'dwnldb'}, true);
        bttn.setRight(pos);
        var handler = function(e) {
            e.stopPropagation();
            e.preventDefault();
            this.pingClick('2-0-Dwnld');
            var cfgIdx = e.target.parentNode.cfgIdx;
            if (cfgIdx == undefined) return;
            this.downloadTblFile([this.m_App.m_Config.TrackConfig[cfgIdx]]);
        }
        bttn = bttn || Ext.get(node.id + 'dwnldb');
        bttn.setListeners({touchstart: handler, click: handler, scope: this});
        return 15;
    },

    addInfoButton: function(gTrk, pos) {
        var node = gTrk.tnEl,
            bId = node.id + 'Info';
            clntLeft = this.bodyDiv.dom.getBoundingClientRect().left,
            tpl = new Ext.Template('<button id="{id}" class="xsv-trackinfo xsv-btntrack" title="Track Info"></button>');
        var bttn = tpl.append(node, {id: bId}, true);
        bttn.setRight(pos);
        var handler = function(e) {
            if (e.button !== 0) return;
            e.stopPropagation();
            e.preventDefault();
            this.pingClick('2-0-Info');
            var cfgIdx = e.target.parentNode.cfgIdx;
            if (cfgIdx == undefined) return;
            this.removeFloatingSelection();
            e.type = 'pointerup'; // to invoke a tooltip without a delay
            this.m_Selection = new SeqView.Selection(this, gTrk.div.dom, [gTrk.area], e, this.bodyDiv.dom.getBoundingClientRect().left);
        }
        bttn = bttn || Ext.get(bId);
        bttn.setListeners({touchstart: handler, click: handler, scope: this});
        return 15;
    },

    addCGVButton: function(node, pos, lnk) {
        if (!lnk.link) return 0; 
        var bId = node.id + 'cgv',
            tpl = new Ext.Template('<button id="{id}" class="xsv-cgvtrack xsv-btntrack" title="Open alignment in Comparative Genome Viewer (CGV)"></button>');
        var bttn = tpl.append(node, {id: bId}, true);
        bttn.setRight(pos);
        var handler = function(e) {
            e.stopPropagation();
            e.preventDefault();
            this.pingClick('2-0-cgv');
            var cfgIdx = e.target.parentNode.cfgIdx;
            if (cfgIdx == undefined) return;
            window.open(lnk.link + this.m_VisFromSeq + '-' + (this.m_VisFromSeq + this.m_VisLenSeq));
        }
        bttn = bttn || Ext.get(bId);
        bttn.setListeners({touchstart: handler, click: handler, scope: this});
        return 15;
    },


    addTrackButtons: function(gTrk, navi, track) {
        var tpl,
            node = gTrk.tnEl.dom,
            pos = 0,
            scrCenter = this.getBodyWidth() >> 1;
        if (navi) { // to avoid overlapping jump buttons
            tpl = new Ext.Template ('<div style="width:' + (scrCenter - 26)
                + 'px;overflow: hidden; white-space: nowrap; text-overflow: ellipsis;">' + node.innerHTML + '</div>');
            node.textContent = '';
            tpl.append(node, {}, true);
        }
        if (!this.m_App.m_NoConfDlg) {
            pos += this.addXButton(node);
            if (track && track.hasOpts) pos += this.addCfgButton(node);
        }
        if (track) {
            if (track.legend_text) pos += this.addInfoButton(gTrk, pos);
            if (!this.m_App.m_noTrackData && SeqView.TM.isDownloadable(track)) pos += this.addDownloadButton(node, pos);
            if (gTrk.area.links) pos += this.addCGVButton(node, pos, gTrk.area.links.find(function(l){return l.label == 'CGV';}) || {});
        }
        if (!navi) return;
        var icon = 'url('+ SeqView.base_url + 'images/jumps.png) ';
        tpl = new Ext.Template('<button class="xsv-jumpbutton" style="left:{left}px;" title="{title}"></button>');
        var lBtn = tpl.append(node, {left: (scrCenter - 26), title: 'Jump left'}, true).dom;
        var rBtn = tpl.append(node, {left: (scrCenter + 8), title: 'Jump right'}, true).dom;
        lBtn.onmouseover = rBtn.onmouseover = function(e) {
            if (this.tmp_gotoFeat || this.m_Loading) return;
            e.target.style.opacity = 0.9;
        }.bind(this);
        lBtn.onmouseout = rBtn.onmouseout = function(e) {
            if (this.tmp_gotoFeat) return;
            var trg = e.target;
            trg.style.width = '18px';
            trg.style.opacity = 0.5;
            trg.style.background = icon + (scrCenter < trg.offsetLeft ? '-18px' : '');
        }.bind(this);
        lBtn.onmouseout({target: lBtn});
        rBtn.onmouseout({target: rBtn});
        lBtn.up = !(rBtn.up = !this.m_Flip);
        lBtn.ontouchstart = rBtn.ontouchstart = lBtn.onclick = rBtn.onclick = function(e) {
            e.stopPropagation();
            e.preventDefault();
            this.pingClick('2-0-3');
            if (this.tmp_gotoFeat || this.m_Loading) return;
            var trg = e.target;
            trg.style.background = SeqView.ExtIconLoading + ' 0px 0px';
            trg.style.width = '11px';
            trg.style.backgroundSize = 'cover';
            this.gotoFeature(this.m_App.m_Config.TrackConfig[trg.parentNode.cfgIdx],
                (trg.up ? 'next' : 'prev') + (e.ctrlKey || e.metaKey ? '_group' : ''),
                function(){ e.target.onmouseout(e); }.bind(this));
        }.bind(this);
    },

    addErrorMsg: function(t, msg) {
        if (t.div != t.tnEl && t.tnEl) {
            var gTrx = this.m_graphTracks
            Ext.removeNode(t.div);
            Ext.removeNode(t.lgEl);
            Ext.each(t.labels, function() { Ext.removeNode(this); });
            var first = t.tnEl.getFirstChild();
            while (t.tnEl.getLastChild() != first) { t.tnEl.removeChild(t.tnEl.getLastChild()); };
            t.div = t.tnEl;
            t.div.applyStyles('top: ' + t.divT + 'px; width: ' + this.getBodyWidth()
                + 'px; height: 18px; position:absolute; font-family:Monospace,Courier New,serif; font-size: 12px;'
                + 'overflow:hidden;white-space: nowrap; text-overflow: ellipsis;'
                + 'background:white;background:linear-gradient(#eeeeee, white, white);');
            t.imgH = t.div.getHeight();
            delete t.area;
            for (var j = idx + 1; j < gTrx.length; j++) {
                var trk = gTrx[j];
                if (!trk || !trk.div) break;
                this.shiftGraphTrack(trk, gTrx[j - 1].div.getBottom(true) - trk.divT);
            }
        }
        delete t.image;
        t.div.dom.cfgIdx = t.cfgIdx;
        t.div.setHtml(t.display_name + ', <a style="color:red;font-size:10px;">' + msg + '</a>');
        this.addTrackButtons(t.div);
    },


    addFailedTrack: function(idx, msg, jobStatus) {
        jobStatus = jobStatus || 'Error';
        var gTrx = this.m_graphTracks,
            tt = [gTrx[idx]];
        if (this.noParallel) {
            tt = this.rqParams.tracks[idx].split('id:');
            tt.shift();
        }
        tt.forEach(function(t) {
            if (typeof t == 'string') {
                t = t.substr(0, t.indexOf(','));
                t = gTrx.find(function(gt) { return gt.id == t; });
                if (!t) return;
            }
            this.addErrorMsg(t,  jobStatus + ': ' + msg);
            if (this.m_App.m_trackEvent) this.m_App.fireEvent('track_image_loaded', this, t.id);
        }, this);
        if (--this.m_Loading == 0) this.applyResultsOnLoad();
        this.checkMolPlacement();
    },

    addTrackImage: function(image, reinitTrack) {
        if (!this.m_Loading) return;
        var gTrx = this.m_graphTracks;
        if (image.idx < 0) {
            var tmpHeight = this.m_topRulerHeight;
            this.m_bottomRulerHeight = this.m_topRulerHeight = image.areas[0].bounds.b + 1;
            var rulTrk = {image: image, divT: 0, imgH: this.m_topRulerHeight, imgT: 0};
            var id = 'trk_Ruler' + this.m_Idx;
            Ext.removeNode(Ext.get(id));
            var topRuler = this.createTrackDiv(rulTrk, id, this.m_vLine).addCls('ruler');
            if (this.m_bottomRulerHeight > 0) {
                rulTrk.imgH = this.m_bottomRulerHeight;
                this.bottomRuler = this.createTrackDiv(rulTrk, id + 'B', topRuler).addCls('ruler');
            }
            if (tmpHeight != this.m_topRulerHeight) this.updateViewHeight(image.idx);

            return this.m_Loading--; //no track image, ruler only
        }
        var szTrackFont = 12,
            tIDs = this.rqParams.tIDs,
            graphDiv = this.bodyDiv,
            trackConfig = this.m_App.m_Config.TrackConfig,
            scrWidth = this.getBodyWidth();

        for (var i in image.areas) { if (!NCBIGBUtils.isNumeric(i)) continue;
            var area = image.areas[i];
            if (!(area.type & (SeqView.AreaFlags.Track | SeqView.AreaFlags.Sequence))) continue;
            
            var track = gTrx[image.idx];
            var cfgTrk = trackConfig[track.cfgIdx] || trackConfig.find(function(t) { return track.id == t.id; }) || {};
            track.cfgIdx = cfgTrk.cfgIdx;

            this.cleanGraphTrack(track);
            if (cfgTrk.legend) {
                Ext.each(cfgTrk.legend, function(){ delete this.bounds; });
                track.legends = cfgTrk.legend;
            }

            Ext.apply(track, {l: area.bounds.l, imgT: area.bounds.t, imgH: image.h - area.bounds.t});
            if (area.type & SeqView.AreaFlags.Sequence) {
                var prms = [area.bounds.t, area.bounds.b];
                //skip this area since we have one more
                image.areas.splice(i, 1);
                area = image.areas[i];
                track.imgT = area.bounds.t;
                track.imgH = image.h - track.imgT;
                this.addSeqPrimers(track, prms);
                this.m_seqPrimers = true;
            }
            track.signature = area.signature.split(';')[0];
            track.display_name = (cfgTrk.key != 'sequence_track' || cfgTrk.check_boxes[0].value)? cfgTrk.display_name : '';
            track.image = image;
    
            if (area.status && area.status.code != 0) track.status = area.status;
            track.area = area;
            area.isSlim = (area.bounds.t == area.bounds.b && track.display_name.length > 0); // slim track (w/o title bar)

            if (track.display_name) { // track name div creation
                var tnEl = track.tnEl = new Ext.Element(document.createElement('div'));
                var displayNameMod = "";
                Ext.each(trackConfig[track.cfgIdx].choice_list, function(ch_list) {
                    if (ch_list.name == "scale") {
                        if (ch_list.curr_value != "linear") {
                            Ext.each(ch_list.values, function(ch_value) {
                                if (ch_value.name == ch_list.curr_value) {
                                    displayNameMod = ch_value.display_name;
                                    return false;
                                }
                            });
                        }
                        return false;
                    }
                });
                Ext.apply(tnEl, {name: track.display_name, parentTrack: track});
                if (displayNameMod && displayNameMod.length > 0) {
                    tnEl.name += " - " + displayNameMod.toLowerCase() + " scaled";
                }
                if (track.status) {
                    tnEl.name += '<div style="color: ' + (track.status.severity == 2 ? 'red' : 'darkgreen')
                               + '; position: relative; top: -5px; font-size:'
                               + (szTrackFont - 2) + 'px; font-family: Arial, san-serif;">' 
                               + (track.status.severity == 2 ? 'Error' : 'Warning')
                               + ': ' + track.status.msg + '</div>';
                }
                tnEl.dom.cfgIdx = track.cfgIdx;
                tnEl.addCls('sv-drag sv-highlight');
                graphDiv.appendChild(tnEl);

                tnEl.applyStyles('text-align: left; position: absolute; top:' + track.divT
                                + 'px; left: 0px; font-family: "Monospace", Courier New, serif; font-size: ' + szTrackFont + 'px;'
                                + (track.area.isSlim ? 'opacity: 0.3; cursor: pointer;' :  ('width:' + scrWidth +'px;'))
                                + 'overflow: hidden; white-space: nowrap; text-overflow: ellipsis;');
                tnEl.insertHtml('afterBegin', tnEl.name);
                if (track.area.isSlim) {
                    track.area.bounds.b += tnEl.getHeight();
                    track.area.tname = tnEl;
                    track.area.bounds.r = tnEl.getWidth();
                    tnEl.dom.area = track.area;
                    tnEl.on('mouseenter',function(e) {
                        this.removeFloatingSelection();
                        if (e.target.area) {
                            if (this.m_tnameDiv) try { this.m_tnameDiv.setStyle('opacity', 0.3); } catch(e) { }
                            this.m_tnameDiv = e.target.area.tname;
                            this.m_Selection = new SeqView.Selection(this, e.target, [e.target.area], e);
                            this.m_tnameDiv.setStyle('opacity', 1);
                        }
                    }, this);
                } else {
                    this.addTrackButtons(track, !(area.type & SeqView.AreaFlags.NoNavigation), cfgTrk);
                    tnEl.on('mouseover',function(e) { this.highlightElement(e, track); }, this);
                }
            }// end of track name div
            if (this.m_App.m_trackEvent) this.m_App.fireEvent('track_image_loaded', this, track.id);
        } //for: image.areas
        // correction for multiple tracks in one image 
        if (track.imgT)
            Ext.each(image.areas, function(a) {
                if (a.type & SeqView.AreaFlags.Ruler) return;
                a.bounds.t -= track.imgT; a.bounds.b -= track.imgT;
            });
         
        for (var i = image.areas.length; i--;) {
            var area = image.areas[i];
            if (!(area.type & SeqView.AreaFlags.Legend)
                || !(area.type & SeqView.AreaFlags.NoSelection)) continue;
            image.areas.splice(i, 1);
            if (!area.id) {
                track.legends = t.legends || [];
                area.opacity = 100;
                track.legends.push(area);
                continue;
            }
            Ext.each(track.legends, function(lg) {
                if (lg.id != area.id) return;
                var c = area.color.split(' ');
                lg.color = '#' + ((1 << 24) + (parseInt(c[0]) << 16) + (parseInt(c[1]) << 8) + parseInt(c[2])).toString(16).slice(1); 
                lg.opacity = parseInt(lg.opacity);
                lg.bounds = area.bounds;
                lg.label = area.signature;
                return false;
            }); 
        }
        // adding comment labels for the picture
        for (i = image.areas.length - 1; i >= 0; i--) {
            if (!(image.areas[i].type & SeqView.AreaFlags.Comment)) continue;
            var area = image.areas.splice(i, 1)[0];
            if (track.id != area.parent_id) continue;
            var newEl = new Ext.Element(document.createElement('div'));
            newEl.insertHtml('afterBegin', area.label);
            var styles = 'white-space:nowrap;position:absolute;font-family:Consoles,Monospace,Courier New, serif; font-size: 11px;';//should always be left:0px;
            if (area.bounds.t) styles += 'top:' + (area.bounds.t + track.divT) + 'px;';
            if (area.bounds.l == -1) styles += 'left:' + area.bounds.r + 'px;';
            if (area.bounds.l == 1) styles += 'right:' + area.bounds.r + 'px;';
            if (area.type & SeqView.AreaFlags.DrawBackground) {
                styles += 'opacity:0.6;background-color:rgb('
                    + (this.m_App.m_Config.Options.curr_color == 'Color' ? '220,220,255' : '195,195,195') + ');';
            }
            if (area.color) styles += 'color:rgba(' + area.color.replace(/ /g, ',') + ');';

            track.labels = track.labels || [];
            track.labels.push(graphDiv.appendChild(newEl.dom));
            newEl.applyStyles(styles);

            if (area.tooltip) var tip = new Ext.ToolTip({target: newEl, html: area.tooltip});
            newEl.on('mouseover', function(e) { this.highlightElement(e, track); }, this);
        }
        track.div = this.createTrackDiv(track, 'trk_' + this.m_Idx + '_' + track.id);
        if (track.legends) track.lgEl = track.div.appendChild(this.createLegendsDiv(track));

        gTrx[image.idx] = track;
        if (!reinitTrack) this.updateViewHeight(image.idx);

        this.m_Loading--;
        this.checkMolPlacement();
    },


    applyResultsOnLoad: function() {
        this.m_App.watchdogStop();
        // log time 
        var ts = new Date().getTime(),
            appTS = this.m_App.m_timeStamps,
            stat = {'SV_graphics_time':  ts - appTS.graphicDataUpdateTS};
        if (appTS.prior2seqconfig) {
            stat['sv-event'] = 'initialization';
            stat['SV_config_time'] = appTS.graphicDataUpdateTS - appTS.prior2seqconfig;
            stat['SV_config_graphics_time'] = ts - appTS.prior2seqconfig;
            delete appTS.prior2seqconfig;
        }
        this.ping(stat);

        var data = this.m_FromCgi,
            tErrors = 0;
            lPos = this.m_Width,
            rPos = 0;

        if (Ext.isIE8) this.m_graphTracks = [];
        var gTrx = this.m_graphTracks;

        // removing tracks hidden during update
        for (var i = gTrx.length; --i >= 0 && gTrx[i].imgH == 0;) this.cleanGraphTrack(gTrx.splice(i, 1)[0]);

        Ext.each(gTrx, function(t, idx) {
            if (!t.image) {
                var name = t.div.getHtml();
                if (name.indexOf('Loading...') != -1) this.addErrorMsg(t, 'Error: track was not loaded');
                return;
            }
            t.div.setLeft(this.m_ScrollPix);
            data.areas = data.areas.concat(t.image.areas || []);
            if (t.image.hairlines && gTrx.length > 1) {
                this.m_featMarkers = this.m_featMarkers || {hls: []};
                t.image.hairlines.sort(function(a, b){return a.pos - b.pos});
                this.m_featMarkers.hls[idx] = {hairlines:  t.image.hairlines, track: t}; 
                lPos = Math.min(lPos, t.image.hairlines[0].pos);
                rPos = Math.max(rPos, t.image.hairlines[t.image.hairlines.length - 1].pos);
            }
            tErrors += t.status; 
        }, this);
        if (tErrors) this.m_App.fireEvent('tracks_errors', gTrx);

        if (this.m_featMarkers) {
            this.m_featMarkers.left = lPos;
            this.m_featMarkers.width = rPos - lPos;
        }

        document.getElementById(this.m_App.m_DivId).style.cursor = '';
        this.hideGeneModeButton(!this.m_showGeneModeButton);
        if (this.m_App.m_Panorama) {
            this.m_Locator = this.m_Locator || new SeqView.Locator(this, this.m_Color, true);
            this.m_App.m_Panorama.m_Locator = this.m_Locator;
        }

        this.m_View.down('#btnReload').enable();

        this.updateViewHeight();
        this.adjustViewWidth();
        this.m_vLine.setHeight(this.m_Height);

        if (this.m_FromCgi.areas && this.m_FromCgi.areas.length > 0) {
            this.m_App.addCustomFeatureFlags({view: this, areas: this.m_FromCgi.areas, callback: this.highlightFlaggedFeatures});
        }
        this.drawFeatMarkers();

        var options = this.rqParams.options;
        if (options.update_locator !== false) this.m_App.updateLocator(this);
        this.pushToHistory();
        this.m_App.updateMarkersPos(this);

        if (this.m_SelectedRangeSet) {
            this.selectRangeSet(this.m_SelectedRangeSet);
            this.m_SelectedRangeSet = null;
        }

        this.blinkAreas(options);
        if (options.callback) options.callback.call(options.cbScope);

        if (!this.m_App.m_DialogShown) this.m_App.resizeIFrame();
        this.m_App.notifyViewLoaded(this);
        this.m_App.updateMarkersSize(this);
        this.m_Loading = 0;
        if (!this.m_App.m_Config.SeqInfo.requested) this.updateStatusInfo();
    },

    createLegendsDiv: function(track) {
        var legEl = new Ext.Element(document.createElement('div'));
        legEl.addCls('sv-drag');
        var bLen = 30;
        var width = this.getBodyWidth();
        var leftOffset = -this.m_ScrollPix;//(this.m_Width - qview_width) >> 1;
        var view = this;
        var img = track.image,
            minT = img.h,
            maxB = 0;
        Ext.each(track.legends, function(lg) {
            if (!lg.bounds) return true;
            minT = Math.min(minT, lg.bounds.t);
            maxB = Math.max(maxB, lg.bounds.b);
        });
        legEl.applyStyles('position: absolute; left:' + leftOffset +'px;'
                        + 'width:' + this.getBodyWidth() +'px; height:' + (maxB - minT) + 'px;'
                        + 'top:' + (track.imgH - maxB + minT - 3) + 'px;');
        Ext.each(track.legends, function(lg, idx) {
            if (!lg.bounds) return true;
            var ib = lg.bounds;
            var tmpl = new Ext.Template('<div class="sv-drag" data-qtip="Adjust subtrack settings" style="left:'
                + ib.l + 'px; top:' + (ib.t - minT)+ 'px;'
                + 'width:' + (ib.r - ib.l) + 'px; height:' + (ib.b - ib.t) +'px; position: absolute;'
                + (lg.parent_id ? '' : 'cursor: pointer') + ';">'
                + '<div class="sv-drag" style="width:' + bLen + 'px; height: 8px; position:relative;'
                + 'background-color:' + lg.color + '; opacity:' + (lg.opacity/100)+ ';"></div>'
                + '<div  class="sv-drag" style="color:blue; font-size: 10px; top: -10px; left:' + (bLen + 5)
                + 'px; width:' + (ib.r - ib.l - bLen - 5) + 'px; overflow: hidden; white-space: nowrap; text-overflow: ellipsis; position:relative;">'
                + (lg.label || lg.signature) + '</div></div>');
            var legTrk = tmpl.append(legEl, {}, true).dom;
            if (!lg.id) return;
            legTrk.onclick = function(e) {
                e.stopPropagation();
                e.preventDefault();
                view.pingClick('4-1', 'Graph_Overlay');
                SeqView.TM.modifyLegendSettings(view, track, idx);
            }
        });
        return legEl.dom;
    },

    createTrackDiv: function(track, id, sibling) {
        var wrapDiv = Ext.get(id);
        if (!wrapDiv) {
            wrapDiv = new Ext.Element(document.createElement('div'));
            wrapDiv.set({id: id});
            if (sibling) wrapDiv.insertAfter(sibling);
            else sibling = this.bodyDiv.insertFirst(wrapDiv, true);
        }
        wrapDiv.applyStyles('top: ' + track.divT + 'px; width: ' + this.m_Width
                          + 'px; height: ' + track.imgH + 'px;position:absolute;overflow:hidden;');

        var trkImg = wrapDiv.getFirstChild();
        if (!trkImg) {
            trkImg = wrapDiv.insertFirst(document.createElement(this.m_trackTag));
        }
        trkImg.set({id: id + this.m_trackTag, class: 'sv-drag sv-highlight sv-dblclick'});
        if (this.m_trackTag == 'img') {
            trkImg.set({src: track.image.img_url,
                style: 'position:relative; overflow:hidden; visibility:visible; draggable:false;'
                + 'top:' + (-track.imgT) + 'px;'});
        } else {
            var cnvs = trkImg.dom;
            cnvs.height = track.imgH;
            cnvs.width = this.m_Width;
            if (track.image)
                cnvs.getContext('2d').drawImage(track.image.img, 0, track.imgT, this.m_Width, track.imgH, 0, 0, this.m_Width, track.imgH);
            else
                NCBIGBUtils.drawRuler(cnvs, track.imgH - 4, 0, this.m_Width, track.start, track.end);
        }
        wrapDiv.dom.style.left = this.m_ScrollPix + 'px';

        return wrapDiv;
    },

    checkMolPlacement: function() {
        var app = this.m_App,
            seqInfo = app.m_Config.SeqInfo;

        if (this.m_Loading > 2 || seqInfo.placements_status != 4 
            || (seqInfo.placements && seqInfo.placements.length)) return;
        var maxRQs = 5;
        seqInfo.placements_status = -4;

        var processError = function(data, text, res) {
            console.log('Molecule placement request error:' + text);
        }

        var processResponse = function(data, text, res) {
            if (data.job_status) {
                if (data.job_status == 'failed' || data.job_status == 'canceled') {
                    processError(null, data.job_status);
                } else {
                    if (maxRQs-- > 0)
                        Ext.defer(SeqView.App.simpleAjaxRequest, 3000, this, [{
                            url: app.m_CGIs.Config + '?job_key=' + data.job_id,
                            context: this,
                            success: processResponse,
                            error: processError}]);
                    else {
                        processError(null, 'Exceeded number of tries');
                        seqInfo.placements_status = 4; // ???EB
                    }
                }
                return;
            }
            if (!data.SeqInfo) return;
            seqInfo.placements = data.SeqInfo.placements;
            seqInfo.placements_status = data.SeqInfo.placements_status; //???EB
            this.updateMolPlacement();
            if (seqInfo.placements && seqInfo.placements.length) this.checkMolPlacement = Ext.emptyFn;
                
        }
        app.AjaxRequest({url: app.m_CGIs.Config + '?req=placements&id=' + app.GI, context: this,
                        success: processResponse,
                        error: processError});

        
    },

    loadFailure: function(res, opts){
        if (res) this.m_App.watchdogReport(res);
        this.m_View.setTitle('Image loading error: ' + res.statusText);
        this.m_Loading--; //??? maybe this.m_Loading = 0
        this.mf_HistoricalBack = false;
    },

    gotoFeature: function(track, dir, callback) {
        callback = callback || function(){};
        this.tmp_gotoFeat = true;
        var vLen = this.m_VisLenSeq,
            vFrom = this.m_VisFromSeq,
            lenSeq = this.m_LenSeq,
            fromSeq = this.m_FromSeq;
        var curMarker = this.m_App.getMarkersInfo().findMarkerByName(SeqView.MarkerNav);
        var sPos = (curMarker) ? curMarker.seq_pos : (vFrom + (vLen >> 1));

        var processError = function(data, text, res) {
            delete this.tmp_gotoFeat;
            console.log('gotoFeature error:' + text);
            callback();
       }

        var processResponse = function(data, text, res) {
            if (data.job_status) {
                if (data.job_status == 'failed' || data.job_status == 'canceled') {
                    processError(null, data.job_status);
                } else {
                    Ext.defer(SeqView.App.simpleAjaxRequest, 500, this, [{
                        url: track.pURL + '?job_key=' + data.job_id,
                        context: this,
                        success: processResponse,
                        error: processError}]);
                }
                return;
            }
            if (typeof data.new_pos !== 'undefined') {
                var dPos = data.new_pos;
                var bOpt = {blinkData: {trackID: track.id, bounds: data.bounds}};
                var newVisFrom = dPos - (vLen >> 1);
                if (newVisFrom < 0) newVisFrom = 0;
                if ((vLen < this.m_App.m_SeqLength) && (dPos < vFrom || dPos > (vFrom + vLen))) {
                    if (newVisFrom < fromSeq || (fromSeq + lenSeq) < (newVisFrom + vLen)) {
                        this.gotoPosRange([dPos, 0], true, bOpt);
                    } else {
                        var shift = Math.round((vFrom - newVisFrom) * this.m_bpPix);
                        if (!this.scrollViewTo(this.m_ScrollPix + shift, false, bOpt))
                            this.blinkAreas(bOpt);
                    }
                } else {
                    this.blinkAreas(bOpt);
                }
                if (curMarker) {
                    curMarker.setSeqPos(dPos, true);
                    curMarker.m_ToolTip.updateMarkerToolTipContent(curMarker.span, dPos);
                }
                else this.m_App.getMarkersInfo().addMarker([dPos], SeqView.MarkerNav, false, 'navy');
            }
            delete this.tmp_gotoFeat;
            callback();
        }
        var data = {track: SeqView.TM.tracksArrayToString(track, true, true),
                   pos: sPos, select: track.selectedSig, dir: dir, id: this.m_App.GI, assm_context: this.m_App.m_AssmContext,
                   navi: 1, from: fromSeq, len: lenSeq, width: this.m_Width};
        if (track.affinity_token) data.affinity_token = track.affinity_token;
        Ext.apply(data, this.m_App.m_Config.visualOptionsUrl());
        this.m_App.AjaxRequest({ url: track.pURL, data: data, context: this,
            success: processResponse,
            error: processError});
    },

    blinkAreas: function(options) {
        if (!options || !options.blinkData) return;
        var area = options.blinkData.area,
            bounds = options.blinkData.bounds,
            tID = options.blinkData.trackID || area.parent_id,
            tpl = new Ext.Template('<div style="position:absolute; border: 1px solid red;"></div>'),
            track;
        if (!tID) return;
        Ext.each(this.m_graphTracks, function(t) { track = t; return t.id != tID; });


        var offset = this.m_FromSeq + (this.m_Flip ? this.m_LenSeq : 0);
        Ext.each(bounds, function(db) {
            db.l = Math.abs(Math.round((db.l - offset) * this.m_bpPix));
            db.r = Math.abs(Math.round((db.r - offset) * this.m_bpPix));
        }, this);
        if (area) bounds = [area.bounds];

        track.blinkData = bounds;
        Ext.each(bounds, function(db) {
            db.elem = SeqView.createSelectionElement(this, track.div.dom, {bounds: db}, tpl);
        }, this);
        function blink(bnum) {
            var stop = bnum-- < 0 || track.blinkData != bounds;
            if (!stop) Ext.defer(blink, 300, this, [bnum]);
            Ext.each(bounds, function(r) {
                if (!stop) r.elem.setVisible(bnum % 2);
                else r.elem.destroy();
            });
        }
        blink(3);
     },

    highlightFlaggedFeatures: function() {
        if (this.m_FromCgi.areas) {
            this.removeFeatMarks();
            Ext.each(this.m_FromCgi.areas, function(a) {
                if( (a.type & SeqView.AreaFlags.Dirty) !== 0) {
                   var sh = new SeqView.SelectionHighlighter(this, a, true);
                   if (!this.m_HighLightedFeature) this.m_HighLightedFeature = [];
                   this.m_HighLightedFeature.push(sh);
                }
            },this);
        }
    },

    updateArea: function(area) {
        var areas = this.m_FromCgi.areas
        for (var i = areas.length; --i >= 0;)
            if (areas[i].signature == area.signature && areas[i].parent_id == area.parent_id) return areas[i];
        return area;
    },

    shiftGraphTrack: function(trk, shift) {
        if (shift == 0) return;
        if (!shift) shift = trk.divT - trk.div.getTop(true);
        else trk.divT += shift;
        trk.div.setTop(trk.divT);
        if (trk.tnEl) trk.tnEl.setTop(trk.divT);
        Ext.each(trk.labels, function() { this.setTop(this.getTop(true) + shift); });
    },

    hideTrack: function(cfgIdx) {
        var gtIdx = this.m_graphTracks.findIndex(function(t) { return t.cfgIdx == cfgIdx; });
        if (gtIdx < 0) return;
        this.m_App.m_actTracks = '';
        var trk = this.m_graphTracks[gtIdx];
        if (trk.image && trk.image.hairlines) this.clearSelectedSig(trk);
        this.removeGraphTracks(gtIdx);
        this.updateViewHeight(gtIdx - 1);

        this.m_App.fireEvent('configuration_changed', this.m_App);

        Ext.each(this.m_FromCgi.areas, function(area, idx) {
            if (area == trk.area || area.parent_id == trk.id) { 
                this.splice(idx, 1);
                return;
            }
        }, this.m_FromCgi.areas, true);

        this.drawFeatMarkers();
        this.m_App.updateMarkersSize(this);
        this.updateStatusBar('tracks', 'Tracks shown: ' + + this.m_graphTracks.length + '/' + this.m_App.m_Config.TrackConfig.length)
        this.m_View.updateLayout();
    },


    shiftDownTrack: function(trkIdx) {
        this.m_App.m_actTracks = '';
        var trx = this.m_graphTracks;
        var src = trx[trkIdx], dst = trx[trkIdx + 1];

        this.shiftGraphTrack(src, dst.imgH);
        this.shiftGraphTrack(dst, -src.imgH);
        trx.splice(trkIdx + 2 , 0, trx[trkIdx]);
        trx.splice(trkIdx, 1);

    },

    clearCache: function() {
        this.m_UrlFrom = this.m_UrlTo = null;
    },

    scrollElements: function(delta) {
        this.m_App.scrollMarkers(this, delta);
        if (this.m_featMarkers) {
            this.m_featMarkers.div.style.left = this.m_featMarkers.left + this.m_ScrollPix + 'px';
        }

        if (this.m_HighLightedFeature) {
            Ext.each(this.m_HighLightedFeature,function(sf) {
                sf.movePix(delta);
            });
        }
        if (this.m_Selection) { this.m_Selection.movePix(delta); }
        var sel_set = this.m_RangeSelectionSet;
        if (sel_set && sel_set.length) {
            var l = sel_set.length;
            for (var x = 0; x < l; ++x)
                sel_set[x].movePix(delta);
        }
    },
    

   updateImagePosition: function(new_pos, drag){
        var delta = this.m_ScrollPix - new_pos;
        if (delta == 0) return false;
        
        delete this.m_ReqLen;
        this.updatingPosition = true;
        var delay = 15;
        var sign = Ext.Number.sign(delta);
        var moveTrx = function(step) {
            step = step || Math.floor((new_pos - this.m_ScrollPix)/10);
            if (step == 0 && new_pos !== this.m_ScrollPix) step = -sign;
            this.m_ScrollPix += step;
            this.scrollElements(-step);
            var left = this.m_ScrollPix + 'px';
            Ext.each(this.bodyDiv.dom.children, function(t){
                if (!t.id.indexOf('trk_')) t.style.left = left;
                
            });
            // returning legends back
            left = -this.m_ScrollPix + 'px';
            Ext.each(this.m_graphTracks, function(t){
                if (t.lgEl) t.lgEl.dom.style.left = left;
            });
            if (this.m_ScrollPix != new_pos) Ext.defer(moveTrx, delay, this);
            else {
                delete this.updatingPosition;
                if (!drag) {
                    this.updateTitle({ from_ui: true, transitional: drag });
                }
                this.m_App.updateLocator(this);
            }
        }
        moveTrx.call(this, drag ? -delta : 0);
        return true;
    },


/*    setTracksOffset: function(offset) {
        offset = offset + 'px';
        Ext.each(this.bodyDiv.dom.children, function(t) {
            if (!t.id.indexOf('trk_')) t.style.left = offset;
        });
    },*/


    scrollViewTo: function(new_pos, drag, options) {
        if (this.m_Loading || this.updatingPosition) return false;

        drag = drag || false;

        this.m_UrlFrom = this.m_UrlTo = null;

        var vis_range = this.toSeq();
        var scrWidth = this.getBodyWidth();
        var bpp = this.m_bpPix;
        var offsetBP = Math.round(new_pos / bpp);
        var basePos = this.m_FromSeq;
        var visFrom = basePos - offsetBP;
        if (this.getFlip()) { 
            basePos += this.m_LenSeq - vis_range[2] - 1;
            visFrom = basePos +  + offsetBP; 
        }
        visFrom = Math.min(Math.max(0, visFrom), this.m_App.m_SeqLength - this.m_VisLenSeq);
        new_pos = Math.round((basePos - visFrom) * bpp);
        if (this.getFlip()) new_pos = Math.floor((visFrom - basePos) * bpp);

        var loading_needed = false;
        if (this.m_Width >= scrWidth) {
            loading_needed = new_pos > 0 || this.m_Width < scrWidth - new_pos;
            if (!this.updateImagePosition(Math.floor(new_pos), drag)) return false;
            this.m_VisFromSeq = visFrom;
            if (!drag) this.pushToHistory(vis_range);
        }
        if (loading_needed) {
            var opt = options || {};
            opt.from_ui = opt.transitional = true;
            this.updateGraphicData(visFrom, this.m_VisLenSeq, opt);
            return true;
        }
        return false;
    },


    updateTitle: function(options) {
        var range = this.toSeq();
        var from = (range[0] + 1);
        var to = (range[1] + 1);

        if (this.getFlip()) {
            var t = from;
            from = to;
            to = t;
        }

        var vp = this.m_App.m_ViewParams;
            title = vp.id + ':&nbsp;' + from.shorten() + '..' + to.shorten() + '&nbsp;(' + range[2].commify() + ' ' + vp.units + ')'; 

        if (this.canFlip() && this.getFlip()) title += ' C';

        var tbar = this.m_View.getDockedItems('toolbar[dock="top"]')[0];
        if (tbar) {
            var tools = tbar.getRefItems();
            // returns true if not all toolbar items are visible 
            var isTBShort = function() { return tbar.el.getRight() < tools[tools.length - 1].el.getRight(); }
            
            var fakeEl = {hide: Ext.emptyFn, show: Ext.emptyFn, setText: Ext.emptyFn};
            var tbSlider = this.m_TbSlider || fakeEl,
                tbSearch = tbar.down('#btnSearch') || fakeEl,
                tbGoTo = tbar.down('#gotoBox') || fakeEl,
                tbFind = tbar.down('#btnFind') || fakeEl,
                tbTools = tbar.down('#btnTools') || fakeEl
                tbDownloads = tbar.down('#btnDownload') || fakeEl,
                tbTracks = tbar.down('#btnTracks') || fakeEl;
            tbSlider.show();
            tbSearch.hide();
            tbGoTo.show();
            tbFind.show();
            tbTools.setText('Tools');
            tbTracks.setText('Tracks');
            tbDownloads.setText('Download');
            var tbTitle = tbar.down('#btnTitle');
            if (tbTitle) {
                var tt = '<b>' + vp.id + ': ' + vp.title
                + '</b><br><br>' + from.commify() + '&nbsp;-&nbsp;' + to.commify() + '<br>'
                + range[2].commify() + '&nbsp;' + vp.units + '&nbsp;shown';
                if (this.canFlip())
                    tt += ',&nbsp;' + (this.getFlip() ? 'reverse complement' : 'positive strand');
                tbTitle.setTooltip(tt);
                tbTitle.setText('<b style="font-size:x-small;">' + vp.id + '</b>');
                if (isTBShort()) tbTitle.setText('');
            }
            //taking care of showing either "Find on Sequence", "Tools", "Tracks" buttons/text
            if (isTBShort()) {
                tbGoTo.hide();
                tbFind.hide();
                tbSearch.show();
            }
            if (isTBShort()) tbTools.setText('');
            if (isTBShort()) tbTracks.setText('');
            if (isTBShort()) tbDownloads.setText('');
            if (isTBShort()) tbSlider.hide('');
        }
        this.m_VisFromSeq = range[0];
        this.m_VisLenSeq = range[2];
        this.updateStatusBar('range', title);

        var flags = {};
        flags.from_ui = options && options.from_ui;
        flags.transitional = options && options.transitional;
        flags.fire = options && options.fire;

        this.m_App.fireEvent('visible_range_changed', this, flags);
    },


    pushToHistory: function(range) {
        if( this.m_InDragAction ) return;
        range = range || this.toSeq(); 
        var flip = this.getFlip();
        var MAX_HISTORY_LENGTH = 30;
        var last = this.m_History[this.m_History.length - 1];

        if (this.mf_HistoricalBack || (last && Math.abs(range[0] - last.from) <= 1 && range[2] == last.len && flip == last.flip) ){
            this.mf_HistoricalBack = false;
            return;
        }
        if (this.m_History.push({from: range[0], len: range[2], flip: flip}) > MAX_HISTORY_LENGTH)
            this.m_History.shift();

        this.mf_HistoricalBack = false;
        this.updateHistoryButtons();
    },

    stepHistory: function(options ){
        this.m_History.pop();
        var prev = this.m_History[this.m_History.length - 1];
        if (prev) {
            this.mf_HistoricalBack = true;
            this.setFlip(prev.flip);
            this.startImageLoading(prev.from, prev.len, options);
        }
        this.updateHistoryButtons();
    },

    updateHistoryButtons: function() {
        var btn_prev = this.m_View.down('#histPrev');
        if (btn_prev) {
            var prev = this.m_History[this.m_History.length - 2];
            btn_prev[prev ? 'enable' : 'disable']();
            var title = 'Back';
            if (prev) {
                var to = prev.from + 1 + (prev.flip ? 0: prev.len);
                var from = prev.from + 1 + (prev.flip ? prev.len : 0);  
                title += ' to ' + from.commify() + '&nbsp;-&nbsp;' + to.commify()
                    + '&nbsp;(' + prev.len.commify() + '&nbsp;' + this.m_App.m_ViewParams.units;
                if (this.canFlip()) title += ',&nbsp;' + (prev.flip ? 'negative' : 'positive');
                title += ")";
            }
            btn_prev.setTooltip(title);
        }
    },

    moveTo: function(vis_from, vis_len, opts){
        if (vis_len == this.m_VisLenSeq) {
            var cur_vis_from_pix = Math.round(this.seq2Pix(this.m_VisFromSeq));
            var new_vis_from_pix = Math.round(this.seq2Pix(vis_from));
            var delta_x = cur_vis_from_pix - new_vis_from_pix;
            if (delta_x == 0) {
                if (opts && opts.callback) opts.callback.call(opts.cbScope);
                return;
            }
            var cur_x = Math.round(this.m_ScrollPix);
            cur_x = cur_x + delta_x;
            var rest_len = this.m_Width + this.m_ScrollPix - this.getBodyWidth();
            if ( cur_x > 0 || (delta_x < 0 && -delta_x > rest_len) ){
                this.startImageLoading(vis_from, vis_len, opts);
            } else {
                // No need to check as it is done above
                this.updateImagePosition(cur_x);
                if (opts && opts.callback) opts.callback.call(opts.cbScope);
                this.updateTitle(opts);
            }
        } else {
            this.startImageLoading(vis_from, vis_len, opts);
        }
    },

    syncTracks: function(tcfg) {
        this.m_graphTracks.forEach(function(t) {
            t.cfgIdx = tcfg.findIndex(function(tc) { return t.id == tc.id; });
            if (t.legends) {
                t.legends.forEach(function(lgTrk) {
                    lgTrk.idx = tcfg.findIndex(function(tc) { return lgTrk.id == tc.id; });
                });
            }
            if (t.tnEl) {
                t.tnEl.dom.cfgIdx = t.cfgIdx;
//                var div = t.tnEl.child('div');
//                if (div) div.setText(t.display_name = tcfg[t.cfgIdx].display_name);
            }
        });
    },

    syncToLocator: function(vis_len) {
        if (!this.m_Locator) return;

        this.moveTo(
            this.m_App.m_Panorama.toSeq(this.m_Locator.getLeft(true)),
            vis_len || this.m_App.m_Panorama.toSeq(this.m_Locator.getWidth()),
            { update_locator: false, from_ui: true }
        );
    },

    checkLocatorWidth: function(width) {
        var len = this.m_App.m_Panorama.toSeq(width+3);
        return SeqView.MinBpp < len/this.getBodyWidth();
    },
    
    getSeqPos: function(x) {
        var leftPos = this.m_VisFromSeq;
        x = Math.floor(x/this.m_bpPix);
        if (this.getFlip()) {
            x = -x;
            leftPos += this.m_VisLenSeq - 1; 
        }
        return leftPos + x;
    },

    hitTest: function(x, y, trk) {
        if (this.m_Loading || this.m_InDragAction)
            return null;

        var the_area = [];
        var isFlip = this.getFlip();
        Ext.each(trk.image.areas, function(area) {
            if (!area.signature) return;
            var bounds = area.bounds;
            var left  = isFlip ? bounds.r : bounds.l;
            var right = isFlip ? bounds.l : bounds.r;
            if (x >= left - 1 && x <= right + 1 && y >= bounds.t && y < bounds.b) {
                the_area.push(area);
//                return false;
            }
        });
//if (the_area.length > 1) { console.log(the_area); }
        return the_area.length ? the_area : null;
    },

    highlightElement: function(e, trk) {
        if (!this.m_App.m_Config.TrackConfig) return;
        var trgt = e.getTarget(this.m_trackTag),
            clientRect = this.bodyDiv.dom.getBoundingClientRect(),
            x = e.getX() - clientRect.left,
            y = e.getY() - Math.round(clientRect.top) - window.pageYOffset;
        if (!trk) { // find current track
            if (!trgt) {
                Ext.each(this.m_graphTracks, function(t, idx) {
                    if (t.divT > y) return false; 
                    trk = t;
                });
                if (!trk) { 
                    console.log('no trk'); return;
                }
            } else {
                if (!this.m_graphTracks.length) return;
                Ext.each(this.m_graphTracks, function(t, idx) {
                    trk = t;
                    return t.div.dom != trgt.parentNode;
                });
            }
        }
        var areas = this.hitTest(x - this.m_ScrollPix, y - trk.divT, trk);
        if (this.m_slimMode && areas && (areas[0].type  & (SeqView.AreaFlags.Track | SeqView.AreaFlags.Sequence))) {
            areas.shift();
            if (!areas.length) return;
            if (this.m_Selection && areas[0].signature != this.m_Selection.areas[0].signature)
                this.removeFloatingSelection();
        }
        if (areas && (!this.m_Selection || this.m_Selection.removed)) {
            this.m_Selection = new SeqView.Selection(this, trk.div.dom, areas, e, clientRect.left);
        } else {
            if (!this.m_Selection || !areas || areas[0].signature != this.m_Selection.areas[0].signature)
                this.removeFloatingSelection();
        }
    },

    onPointerDown: function(e) {
        if (this.m_ContextMenu) this.m_ContextMenu.destroy();
        if ((e.pointerType == 'mouse' && e.button) || this.m_Loading) return;
        this.m_XY = e.getXY();

        if (this.updatePosTT(e)) {
            // Add range selection
            e.stopEvent();
            if (!e.ctrlKey) this.removeSelectionsWithNoPinnedToolTips();
            this.m_InDragAction = true;
            new SeqView.RangeSelection(this, [this.m_XY[0]], e);
            this.m_XY = null;
//            this.bodyDiv.setStyle('cursor', 'e-resize');
            return;
        }
        if (e.pointerType != 'mouse')
            this.m_deferredContext = Ext.defer(this.showContextMenu, 2000, this, [this.m_XY]);

        this.time = (e.pointerType == 'touch') ? e.getTime() + 200 : 0;

        if (Ext.fly(e.getTarget()).hasCls('sv-drag')) {
            this.m_InDragAction = true;
            e.stopEvent();
            this.bodyDiv.setStyle('cursor', 'move');
            this.m_App.eventHandlers('on', this);
        }
    },

    onPointerUp: function(e) {
        this.m_InDragAction = false;
        if (!this.m_XY) return;
        if (this.time > e.getTime()) {
            this.removeFloatingSelection();
            this.highlightElement(e);
        }
        if (this.m_deferredContext) {
            clearTimeout(this.m_deferredContext);
            this.m_deferredContext = 0;
        }
        e.stopEvent();
        this.bodyDiv.setStyle('cursor', '');//'default');

        if (this.m_trackClone) {
            Ext.removeNode(this.m_trackClone.dom);
            delete this.m_trackClone;
            delete this.m_App.m_currentTrackSetId;
            this.m_App.fireEvent('configuration_changed', this.m_App);
            this.drawFeatMarkers();
            this.pingClick('2-0-0');
        }
        if (this.m_scrolled) {
            this.updateTitle({from_ui: true, transitional: false, fire: true });
            delete this.m_scrolled;
            clearTimeout(SeqView.scrolled);
            SeqView.scrolled = setTimeout(function() { delete SeqView.scrolled; }, 500);
            if (this.m_ReqFrom != this.m_VisFromSeq) this.pingClick('2-0-2');
        }
        this.m_XY = null;
        this.m_App.eventHandlers('un', this);

        this.pushToHistory();
    },

    onPointerMove: function(e) {
        NCBIGBUtils.ClearBrowserSelection();
        if (this.m_Loading) return;
        if (this.m_InDragAction) {
            if (!this.m_XY) return;
            if (e.button > 0) {
                this.onPointerUp(e)
                return;
            }

            var sensivity = 2;
            var xy = e.getXY();
            var delta_x = this.m_XY[0] - xy[0];
            var delta_y = this.m_XY[1] - xy[1];
            if (delta_x == 0 && delta_y == 0) return;
            var XminusY = Math.abs(delta_x) - Math.abs(delta_y);
            if (Math.abs(XminusY) <= sensivity) return;
            if (this.m_deferredContext) {
                clearTimeout(this.m_deferredContext);
                this.m_deferredContext = 0;
            }

            if (XminusY > 0 ) {
                var cur_x = this.m_ScrollPix;
                cur_x = cur_x - delta_x;
                this.m_XY = xy; // save new values
                this.m_scrolled = true;
                this.scrollViewTo(cur_x, true);
            } else { // tracks moving processing
                var bodyY = this.bodyDiv.getY(),
                    tracks = this.m_graphTracks,
                    pY = this.m_XY[1] - bodyY,
                    src = -1, dst;
                var tClone = this.m_trackClone;
                if (tClone) {
                    if (!tracks.some(function(t) { src++; return t.id == tClone.srcId; })) src = -1;
                } else {
                    if (!tracks.some(function(t) { src++; return pY > t.divT && pY < t.divT + t.imgH; })) src = -1;
                }
                if (src >= 0) {
                    if (!tClone) {
                        tClone = this.createTrackDiv(tracks[src], 'trk_clone', this.bodyDiv.dom.childNodes[tracks.length]);
                        tClone.dom.style.border = 'green solid 1px';
                        tClone.srcId = tracks[src].id;
                        tClone.mouseOffsetY = pY - tracks[src].divT;
                    }
/*                    if (e.ctrlKey) 
                        console.log('Math.min(Math.max(' + this.m_topRulerHeight + ',' + xy[1] + '-' + bodyY + '-' + tClone.mouseOffsetY +'),' +
                            this.bodyDiv.getBottom(true) + '-' + tClone.getHeight()+'-' + this.m_bottomRulerHeight +'))');*/
                    tClone.setLocalY(Math.min(
                        Math.max(this.m_topRulerHeight, xy[1] - bodyY - tClone.mouseOffsetY),
                        this.bodyDiv.getBottom(true) - tClone.getHeight() - this.m_bottomRulerHeight));
                    var top = tClone.getTop(true),
                        bottom = tClone.getBottom(true);
                    if (delta_y > 0) {// moving up
                        for (dst = src; dst > 0; dst--) {
                        if (e.ctrlKey) console.log(src+'/'+dst + ': ' + tracks[dst].divT + '-' + top , tracks[dst].divT - top);
                        if (tracks[dst].divT - top < 10)  break; }
                    } else {// moving down
                        for (dst = src; dst < tracks.length - 1; dst++) { if (tracks[dst].divT >= top) break; }
                    }
                    if (e.ctrlKey) console.log(dst);
                    if (dst != src) {
                        for (var i = Math.min(src, dst), max = Math.max(src, dst); i < max; i++) {
                            this.shiftDownTrack(i);
                        }
                        this.m_XY = xy;
                        for (var i = 0; i < tracks.length; i++) {
                            this.m_App.m_Config.TrackConfig[tracks[i].cfgIdx].order = i;
                        }
                    } 
                }
                this.m_trackClone = tClone;
            }
        } else {
            if (!this.updatePosTT(e)){//, this.m_App.m_hiddenOptions.indexOf('nopostt') >= 0)) {
                 this.bodyDiv.setStyle('cursor', '');
                 if (e.event.movementY == undefined || Ext.isSafari || (e.event.movementX | e.event.movementY)) {
                    var trgt = Ext.fly(e.getTarget());
                    if (trgt.hasCls('sv-highlight')) this.highlightElement(e, trgt.parentTrack);
                }
            }
        }
        e.stopEvent();
    },

    checkRuler: function(e) {
        return e.getTarget('div.ruler') || !e.target.id.indexOf('line') || (!e.getTarget('div.marker_label') && e.getLocalXY()[1] < this.m_topRulerHeight);
    },

    onContextMenu: function(e) {
        e.stopEvent();
        this.showContextMenu(e.getXY());
    },

    updatePosTT(e) {
        var tt = this.m_posTT;
        if (this.checkRuler(e) || e.ctrlKey) {
            var xy = e.getXY();
            this.m_vLine.setX(xy[0]);
            this.m_vLine.setStyle('display', 'block');
            tt.setHtml(this.m_App.posToLocalDisplay(this.pix2Seq(xy[0] - this.bodyDiv.getX() - this.m_ScrollPix)));
            xy = [xy[0] - 5, xy[1] + (this.m_App.m_hiddenOptions.indexOf('bottomtt') < 0 ? -30 : 20)];
//            tt.pageX = tt.x = xy[0]; tt.pageY = tt.y = xy[1]; 
            if (!tt.isVisible()) {
//            if (!this.checkRuler(e)) console.log(e.target);
                tt.setDisabled(false);
//                tt.setTarget('');
                tt.show();
//                tt.setTarget(this.bodyDiv);
            }
            tt.setXY(xy);
//                          console.log(this.m_posTT.pageX, this.m_posTT.x, xy[0] - 5, xy[1] + (this.m_App.m_hiddenOptions.indexOf('bottomtt') < 0 ? -30 : 20));
            return true;
        }
        tt.setDisabled(true);
        tt.hide();
        this.m_vLine.setStyle('display', 'none');
        return false;        
    },

    // callback allows to customize the menu if you clicked on an object
    showContextMenu: function(xy, feat_pos, callback) {
        if (this.m_Loading) return;
        if (this.m_deferredContext) {
            clearTimeout(this.m_deferredContext);
            this.m_deferredContext = 0;
        }
        var menu = new Ext.menu.Menu();
        var x_pos = xy[0] - this.bodyDiv.getX();
        var seq_pos = this.pix2Seq(x_pos - this.m_ScrollPix);
        menu.add(
            {text:'Select a Range', icon: SeqView.base_url + 'images/range.png', scope:this,
             handler:function() { this.pingClick('2-1-Range'); this.selectRange([seq_pos, this.pix2Seq(x_pos + 50 - this.m_ScrollPix)]); }
        });

        menu.add(
            {text:'Set New Marker At Position', iconCls:'xsv-markers', scope:this,
             handler:function() {this.pingClick('2-1-0'); this.m_App.newMarkerDlg(this, x_pos); }
        });
        if (this.m_RangeSelectionSet  &&  this.m_RangeSelectionSet.length === 1) {
            menu.add(
                {text:'Set New Marker On Selected Range',iconCls:'xsv-markers', scope:this,
                 handler:function() {
                    this.pingClick('2-1-1');
                    this.m_App.addMarker(this.m_RangeSelectionSet[0].range);
                    this.removeRangeSelection();
                }
            });
        }
        menu.add('-');
        if (this.m_App.m_Origin) {
            menu.add({text:'Reset Sequence Origin', iconCls:'xsv-origin', scope:this,
                handler:function() {this.m_App.clearOrigin(); this.pingClick('2-1-2');}});
        } else {
            menu.add({text:'Set Sequence Origin At Position', iconCls:'xsv-origin', scope:this,
                handler:function() { this.m_App.setOrigin(this, x_pos); this.pingClick('2-1-2');}});
        }
        if (feat_pos !== undefined)
            menu.add({text:'Set Sequence Origin At Feature', iconCls:'xsv-origin', scope: this,
                          handler:function() { this.m_App.showOriginDlg(feat_pos);} });


        if (this.canFlip()) {
            menu.add({
                text:'Flip Sequence Strands', iconCls:'xsv-flip-strands', scope: this,
                enableToggle: true, pressed: this.getFlip(),
                handler:function() { this.flipStrand(); }
            });
        }

        menu.add({iconCls:'xsv-zoom_plus', text:'Zoom In', scope: this,
            handler:function() {this.zoomIn(seq_pos); this.pingClick('2-1-3');}});
        menu.add({iconCls:'xsv-zoom_minus', text:'Zoom Out', scope: this,
            handler:function() {this.zoomOut(seq_pos); this.pingClick('2-1-4');}});
        menu.add({iconCls:'xsv-zoom_seq', text:'Zoom To Sequence', scope: this,
            handler:function() {this.zoomSeq(seq_pos); this.pingClick('2-1-5');}});
        menu.add({iconCls:'xsv-zoom_range', text:'Zoom On Range', scope: this,
            handler:function() {this.zoomRange(); this.pingClick('2-1-6');}});
/*        if (this.m_App.mf_MultiPanel) {
            menu.onBefore('show', function(menu) { 
                menu.down('#addRPanel')[this.m_RangeSelectionSet.length ? 'enable' : 'disable']();
            }, this);
            menu.add({iconCls:'xsv-new_view', text:'Add New Panel on Range', itemId: 'addRPanel',
                tooltip:'Create New Graphical Panel on Selected Range',
                scope:this, handler:function() { this.addNewPanel('2-1-7', this.getTotalSelectedRange()); }
            });
        }
*/
        menu.add('-');
        // Tools - Blast and Primer Blast
        var noPBlast = (this.m_App.m_ViewParams['acc_type'] !== 'DNA');
        var tools_submenu = new Ext.menu.Menu();
        if (this.m_RangeSelectionSet  &&  this.m_RangeSelectionSet.length > 0  &&  this.m_RangeSelectionSet.length <= 2) {
            tools_submenu.add({iconCls:'xsv-blast', text:'BLAST Search (Selection)', scope: this,
                handler:function() { this.blastSelection(); this.pingClick('2-1-8-1');} });
            tools_submenu.add({iconCls:'xsv-primer', text:'Primer BLAST (Selection)', scope: this, disabled: noPBlast,
                handler:function() { this.primerBlast(); this.pingClick('2-1-8-2');} });
        }
        tools_submenu.add({text:'BLAST Search (Visible Range)', iconCls:'xsv-blast', scope: this, handler:function() {
                              this.m_App.blast(this, [this.m_VisFromSeq, this.m_VisFromSeq + this.m_VisLenSeq - 1]);
                              this.pingClick('2-1-8-3');
                        } });
        tools_submenu.add({text:'Primer BLAST (Visible Range)', iconCls:'xsv-primer', scope: this, disabled: noPBlast, handler:function() {
                              this.primerBlast(false, [ [this.m_VisFromSeq, this.m_VisFromSeq + this.m_VisLenSeq - 1] ]);
                              this.pingClick('2-1-8-4');
                        } });
        menu.add({text:'BLAST and Primer Search', tooltip:'Tools', iconCls:'xsv-views_tools', menu:tools_submenu});

        menu.add(
            {text:'Printer-Friendly (PDF/SVG)',iconCls:'xsv-printer', disabled: this.m_App.m_NoPDF, scope: this,
             handler: this.downloadImgFile}
        );
        if (this.m_App.m_NoConfDlg !== true) {
            menu.add('-', {iconCls: 'xsv-config', text: 'Configure tracks',
                hidden: this.m_App.countActiveSeqPanels() != 1,
                scope: this, handler:function() { this.m_App.showTracksConfigDlg('2-1-cfg'); }
            });
        }

        if (callback) callback(menu, download_submenu);

        if( menu.items.length > 0 ){
            menu.showAt(xy);
            this.m_ContextMenu = menu;
        }
    },


    downloadData: function(isRangeSelected, format, range) {
        if (isRangeSelected && this.m_RangeSelectionSet.length == 0) {
            Ext.Msg.show({ title:'Download', msg: 'Error: no selection is found', buttons: Ext.Msg.OK, icon: Ext.Msg.ERROR });
            return;
        }
        var ranges = "";
        var first = this.m_App.m_Config.SeqInfo.length, last = 0;
        var flip = this.m_App.m_Flip;
        var strand = flip ? "(-)" : "";

        if (range && range.length > 0) {
            if (!flip)
                ranges = range[0] + "-" + range[1];
            else
                ranges = range[1] + "-" + range[0];
            first = range[0];
            last = range[1];
        }
        else if (!range && isRangeSelected) {
            for (x = 0; x < this.m_RangeSelectionSet.length; x++) {
                range = this.m_RangeSelectionSet[x].range;
                if (x > 0 ) ranges += ",";
                if (!flip)
                    ranges += range[0] + "-" + range[1];
                else
                    ranges += range[1] + "-" + range[0];
                last = Math.max(range[1], last);
                first = Math.min(range[0], first);
            }
        } else {
            last = this.m_VisFromSeq + this.m_VisLenSeq - 1;
            first = this.m_VisFromSeq;
            if (!flip)
                ranges = first + "-" + last;
            else
                ranges = last + "-" + first;
        }
        var params = '?id='+ this.m_App.GI + '&netcache=true&format=' + format + '&ranges=' + ranges,
            fname = '&filename=' + this.m_App.m_Config.SeqInfo.id + '[' + (++first) + '..' + (++last) + ']' + strand + '.' + (format != 'flat' ? 'fa' : format);

        if (ranges.length > 0) this.m_App.requestFile(this, params, fname);     
    },


    downloadAllMarkers: function (format) {
        var mi = this.m_App.getMarkersInfo();
        if (!mi.hasMarkers()) {
            Ext.Msg.show({ title: 'Download', msg: 'Error: no markers found', buttons: Ext.Msg.OK, icon: Ext.Msg.ERROR });
            return;
        }
        var ranges = '';
        var first = this.m_App.m_Config.SeqInfo.length, last = 0;
        var flip = this.m_App.m_Flip;

        this.m_App.m_MarkersInfo.forEachMarker(function(m) {
            var end = m.seq_pos + (m.span ? m.span - 1 : 0);
            if (ranges) ranges += ',';
            if (m.marker_name && m.marker_name != '')
                ranges += m.marker_name + ":";
            if (!flip)
                ranges += m.seq_pos + "-" + end;
            else
                ranges += end + "-" + m.seq_pos;
            last = Math.max(end, last);
            first = Math.min(m.seq_pos, first);
        });
        var params = '?id=' + this.m_App.GI + '&format=' + format + '&ranges=' + ranges,
            fname = '&filename=' + this.m_App.m_Config.SeqInfo.id + '[' + (++first) + '..' + (++last) + ']' + (flip ? '(-).' : '.') + (format != 'flat' ? 'fa' : format);

        if (ranges.length > 0) this.m_App.requestFile(this, params, fname);     
    },

    downloadTblFile: function(tracks) {
        var app = this.m_App,
            gview = this,
            formPanel,
            includeRNAOverrule = true,
            maxSupportedIncludeRNA = 1000000000;
        var tracks_list = Ext.create('Ext.grid.Panel', {
            store: tracks,
            columns: [{ text: '<b>Select one track</b> <i>(choose from following supported tracks in view)</i>', dataIndex: 'display_name', hideable: false, flex: 1 }],
            height: 150, width: 465,
            hidden: tracks.length == 1,
            renderTo: Ext.getBody(),
            listeners: { scope: this,
                select: function(sm, rec) {
                    formPanel.initForm(rec.data, this.m_tableFormat);
                }
            }
        });

        var reportStatus = function(text, icon) {
            var st = win.down('#status');
            if (typeof text == 'string') st.setText(text);
            st.setStyle('background', icon || '');
            st.setStyle('padding', (st.el.getStyle('background-image') != 'none' ? '0 22px' : '0 3px'));;
        };

        var updateLinkToSource = function(track)
        {
            var LTS_ed =  formPanel.down('#linkToSource');
            LTS_ed.setHidden(true);
            if (track.data_provider.indexOf('_rmt' > 0))
            {
                var url = track.annots.find(function(a) { return a.search(/(https?|ftp):\/\//) == 0; });
                if (url )
                {
                    LTS_ed.setText(' <br>Note: The downloaded VCF will only include the data used to render this track and corresponding tooltips. You can download the original VCF in its entirety here: <a href='
                    + url + ' target=\"_blank\">' + url.split('/').pop() + '</a><br>',false);
                    LTS_ed.setHidden(false);
                }
                LTS_ed._rmt = true;
                var tf = formPanel.down('#tdType').getValue().rb || 'GFF3';
                LTS_ed.setHidden(tf!='VCF')

            }
            else
                LTS_ed._rmt = false;
        }

        var updateIncludeRNACheckbox = function(rWidth) {
            var RNA_cb =  formPanel.down('#includeRNA'),
                LTS_ed =  formPanel.down('#linkToSource'),
                showLinkToSrc = true,
                boxLabel = 'Include RNA and CDS features ';
            this.m_includeRNA = RNA_cb.checked;
            this.m_tableFormat = formPanel.down('#tdType').getValue().rb || 'GFF3';
            switch (this.m_tableFormat) {
                case 'GFF3':
                    includeRNAOverrule = (rWidth > maxSupportedIncludeRNA);
                    if (!includeRNAOverrule) {
                        if (RNA_cb.disabled) {
                            RNA_cb.setDisabled(false);
                            RNA_cb.setValue(this.m_includeRNALastKnownValue);
                        }
                    } else {
                        if (!RNA_cb.disabled) this.m_includeRNALastKnownValue = this.m_includeRNA;
                        RNA_cb.setDisabled(true);
                        RNA_cb.setValue(false);
                        boxLabel +='(not supported for regions larger than 10Mbp)';
                    }
                    showLinkToSrc = false;
                    break;
                case 'CSV':
                    includeRNAOverrule = true;
                    // break is missed intentionaly
                case 'BED':
                    if (!RNA_cb.disabled) this.m_includeRNALastKnownValue = this.m_includeRNA;
                    RNA_cb.setDisabled(true);
                    RNA_cb.setValue(false);
                    boxLabel += '(currently not supported for this format)'
                    showLinkToSrc = false;
                    break;
           }
           RNA_cb.setBoxLabel(boxLabel);

           if (LTS_ed._rmt)
           {
                LTS_ed.setHidden(!showLinkToSrc)
           }
        }.bind(this);

        this.m_includeRNA = this.m_includeRNA !== false;
        this.m_includeRNALastKnownValue = true;

        var downloadData = function() {
            this.m_App.handlePos(formPanel.down('#range').getValue(), { allow_equal: true, ask_user: true, scope: this,
                failure: function(err_msg, options) { Ext.MessageBox.alert('Error', err_msg); },
                success: function(pos_range, options) {
                    if (pos_range.length !== 2) {
                       options.failure.call(options.scope, "Range required", options);
                       return;
                    }
                    formPanel.down('#save').disable();

                    updateIncludeRNACheckbox(pos_range[1] - pos_range[0]);
                    reportStatus('File creation started', SeqView.ExtIconLoading + ' no-repeat');
                    this.pingClick('12-1-' + this.m_tableFormat);
                    this.loader = SeqView.App.downloadTracks({
                        id: this.m_App.m_ViewParams.id,      // sequence accession
                        tracks: [tracks_list.getSelection()[0].data],    // tracks as objects
                        ranges: [pos_range], // 0-based
                        //wait: 10000,         // debug: wait 30 sec before pass data
                        format: this.m_tableFormat,
                        assm_context:this.m_App.m_AssmContext,
                        transcr_prot_features: includeRNAOverrule ? false : this.m_includeRNA
                        }, Ext.emptyFn);
                    this.loader.then(function () {
                        reportStatus('File created');
                        delete gview.loader;
                        formPanel.down('#save').enable();
                    }).catch(function (e) {
                        Ext.MessageBox.alert('Error', e.message || 'Unknown error occured');
                        delete gview.loader;
                        reportStatus('No File created');
                        formPanel.down('#save').enable();
                    });
                 }
            });
        }.bind(this);

        var btnHandler = function(btn, event) {
            switch (btn.itemId) {
                case 'help':
                    gview.pingClick('12-4');
                    window.open('/tools/sviewer/seqtrackdata/');// + '&inline=true');
                    break;
                case 'save': downloadData();
                    break;
                case 'cancel':
                    if (gview.loader != undefined) {
                        gview.loader.cancel();   // call loader cancel
                        delete gview.loader;
                        gview.pingClick('12-3');
                    } else gview.pingClick('12-2');
                    win.close();
                    break;
            }
        }

        var win = new Ext.Window({
            title: tracks_list.hidden ? tracks[0].display_name : 'Download Track Data',
            constrain: true,
            app: this.m_App,
            modal: true,
            layout: 'fit',
            minWidth: 500, width: 500,
            minHeight: tracks_list.hidden ? 310 : 510,
            //height: 'auto',//510,
            plain: true,
            cls: 'SeqViewerApp',
            bbar: [{ iconCls: 'x-tbar-loading', xtype: 'tbtext', height:16, itemId: 'status', html: 'Ready' }]
        });
        var adjustedRange1 = this.m_App.posToLocalDisplay(this.m_VisFromSeq);
        var adjustedRange2 = this.m_App.posToLocalDisplay(this.m_VisLenSeq + this.m_VisFromSeq - 1);
        var rangeWidth = adjustedRange2 - adjustedRange1;
        includeRNAOverrule = (rangeWidth > maxSupportedIncludeRNA);
        formPanel = new Ext.FormPanel({
            labelWidth: 1, // label settings here cascade unless overridden
            frame: true,
            bodyStyle: 'padding:5px 5px 0',
            items: [
                { xtype:'fieldset', title: 'Download range', autoHeight: true, defaultType: 'textfield',
                  items :[
                        { xtype: 'label', html: 'Possible range formats include 10K-20K, 10:20, 20000-30000, 5 to 515, 1...246<br><br>'},
                        { itemId: 'range', width: '100%', value: adjustedRange1 + ':' + adjustedRange2,
                          listeners: {
                            scope: this,
                            focusleave: function (lb, nv, ov) {
                                this.m_App.handlePos(formPanel.down('#range').getValue(), {
                                    allow_equal: true,
                                    ask_user: false,
                                    success: function (pos_range, options) {
                                        if (pos_range.length !== 2) return;
                                        adjustedRange1 = pos_range[0];
                                        adjustedRange2 = pos_range[1];
                                        rangeWidth = adjustedRange2 - adjustedRange1;
                                        updateIncludeRNACheckbox(rangeWidth);
                                    },
                                    failure: function (err_msg, options) { Ext.MessageBox.alert('Error', err_msg); }
                                });
                            }
                        }}]
                },
                tracks_list,
                { xtype: 'displayfield', itemId: 'formats', value: 'Choose a format:',
                  style: {color:'grey', "text-align": 'left', "margin-top":'10px', "margin-left":'6px'}},,
                { xtype: 'radiogroup', columns: 1, vertical: true, itemId: 'tdType',
                  defaults: {
                    name: 'rb',
                    listeners: { change: function(rf, val) { if (val) updateIncludeRNACheckbox(adjustedRange2-adjustedRange1); } } 
                  }
                },
                { xtype:'checkbox', itemId: 'includeRNA', checked: includeRNAOverrule ? false : this.m_includeRNA,
                  boxLabel: 'Include RNA and CDS features' + (includeRNAOverrule ? '(Not supported for regions larger than 10Mbp)' : ''),
                  disabled: includeRNAOverrule }
                  ,
                  { xtype:'label', itemId: 'linkToSource', hidden: true  }
            ],


            buttons: [
                { text:"Help", tooltip:'Help',ariaLabel:'Help', itemId: 'help', dock:'left', scope: this, disabled: false, handler: btnHandler },
                '->',
                { text: 'Download', itemId: 'save', scope: this, disabled: false, handler: btnHandler },
                { text: 'Cancel',itemId: 'cancel', scope: this, handler: btnHandler
            }]
        });
        formPanel.initForm = function(track, curFmt) {
            if (this.curTrackKey != track.key) {
                var fmtList = ['GFF3','CSV','BED', 'VCF'];
                var formatCfg = {
                    gene_model_track: {GFF3: true, CSV: true, BED: true, default: 'GFF3'},
                    aggregate_feature_track: {GFF3: false, CSV: false, BED: true, default: 'BED'},
                    feature_track: {GFF3: false, CSV: false, BED: true, default: 'BED'},
                    SNP_track: {GFF3: false, CSV: false, BED: true, VCF: true, default: 'BED'},
                    dbvar_track: {GFF3: false, CSV: false, BED: true, VCF: true, default: 'BED'}};
                //if (!curFmt || !formatCfg[track.key][curFmt]) 
                curFmt = formatCfg[track.key].default;
                
                var fmts = this.down('#tdType');
                fmts.removeAll();
                fmtList.forEach(function(fmt) {
                    fmts.add({boxLabel: fmt, inputValue: fmt, disabled: !formatCfg[track.key][fmt]});
                });
                fmts.setValue({rb: curFmt});
            }
            this.curTrackKey == track.key;
            updateLinkToSource(track);
        }
        win.add(formPanel);
        tracks_list.getSelectionModel().select(0);
        win.show();
        updateIncludeRNACheckbox(rangeWidth);


    },

    saveFile: function(fileURL, pingarea) {
        if (pingarea) this.pingClick(pingarea);
        var form = Ext.DomHelper.append(document.body, { tag: 'form', method: 'post', action: fileURL});
        document.body.appendChild(form);
        form.submit();
        document.body.removeChild(form);
    },

    downloadImgFile: function() {
        var formPanel,
            app = this.m_App;
        var btnHandler = function(btn, event) {
            app.showMessage('', true);
            if (!formPanel.fileURL) {
                this.pingClick('11-' + (this.m_App.m_PDFDownloadType ? 1 : 2));
                app.handlePos(formPanel.down('#range').getValue(),
                    { nextAction: btn.itemId, allow_equal: false, ask_user: true, formPanel: formPanel, scope: this,
                      target: this.m_App.m_PDFDownloadType ? 'pdf' : (this.m_App.m_hiddenOptions.indexOf('svgz') >= 0 ? 'svgz' : 'svg'), 
                      success: this.createImageFile,
                      failure: function (error, options) { app.showMessage(error, true); }}
                );
                return;
            }
            switch (btn.itemId) {
                case 'save':
                    this.saveFile(formPanel.fileURL, '11-3');
                    break;
                case 'view':
                    window.open(formPanel.fileURL + '&inline=true');
                    break;
            }
        }

        var win = new Ext.Window({
            title: 'Download Image',
            constrain: true,
            app: this.m_App,
            modal: true,
            layout:'fit',
            minWidth:500,
            width:500,
            height:330,
            plain: true,
            cls: 'SeqViewerApp',
            bbar: [{ xtype: 'tbtext', itemId: 'status', id: 'status' + this.m_App.m_Idx, html: 'Ready' }]
        });
        var adjustedRange1 = this.m_App.posToLocalDisplay(this.m_VisFromSeq);
        var adjustedRange2 = this.m_App.posToLocalDisplay(this.m_VisLenSeq + this.m_VisFromSeq - 1);
        var selfSVA = this.m_App;

        if(this.m_App.m_PDFDownloadType===undefined)
            this.m_App.m_PDFDownloadType = true;
        var fpHandler = function() {
            formPanel.down('#save').show();
            delete formPanel.fileURL;
        }

        formPanel = new Ext.FormPanel({
            labelWidth: 1, // label settings here cascade unless overridden
            frame:true,
            bodyStyle:'padding:5px 5px 0',
            items: [
            { xtype:'fieldset', title: 'Enter Sequence Range', autoHeight: true, defaultType: 'textfield',
                items :[
                    { xtype: 'label', html: 'Possible range formats include 10K-20K, 10:20, 20000-30000, 5 to 515, 1...246<br><br>' },
                    { itemId: 'range', width: '100%', value: adjustedRange1 + ':' + adjustedRange2, listeners: { change: fpHandler }}
                ]},
            { xtype:'checkbox', itemId: 'title', checked: (this.m_App.m_PDFtitle) ? this.m_App.m_PDFtitle : false, boxLabel: 'Include Title',
                listeners: { change: fpHandler }},
            { xtype:'checkbox', itemId: 'simplified', checked: (this.m_App.m_PDFcomp) ? this.m_App.m_PDFcomp : false, boxLabel: 'Simplified color shading (allows greater compatibility with image editors)',
                listeners: { change: fpHandler }},
            { xtype:'displayfield', id: 'sv-uplerr' + this.m_App.m_Idx, value: '', hidden: true,
                style: {color:'red', "text-align": 'left', "margin-top":'10px', "margin-left":'6px'}},
            { xtype: 'displayfield', id: 'sv-vector' + this.m_App.m_Idx, value: 'File type: (PDF and SVG contain vector graphics for high quality images)',
                 style: {color:'grey', "text-align": 'left', "margin-top":'10px', "margin-left":'6px'}},
            { xtype: 'radiogroup', columns: 1, vertical: true,
              listeners: {
                  change: function (cb, nv, ov) {
                      selfSVA.m_PDFDownloadType = nv.isPDF;
                      formPanel.down('#view')[nv.isPDF ? 'show' : 'hide']();
                      delete formPanel.fileURL;
               }},
               items: [
                    { boxLabel: 'PDF', name: 'isPDF', inputValue: true, checked: this.m_App.m_PDFDownloadType },
                    { boxLabel: 'SVG', name: 'isPDF', inputValue: false, checked: !this.m_App.m_PDFDownloadType }
                ]}
            ],
            buttons: [{
                text: 'Preview', itemId: 'view', scope: this, disabled: false, hidden: !this.m_App.m_PDFDownloadType, handler: btnHandler },{
                text: 'Download', itemId: 'save', scope: this, disabled: false, handler: btnHandler },{
                text: 'Cancel', scope: this, handler: function() { win.close(); }
            }]
        });
        win.add(formPanel);
        win.show();
    },

    createImageFile: function(pos_range, options) {
        var app = this.m_App,
            width = this.getBodyWidth();
        var reportStatus = function(text){
            var st = Ext.getCmp('status' + options.scope.m_App.m_Idx);
            st.setHtml(text);
        }
        if (width > 8000) {
            var msg = 'Image is too wide (projected  width is ' + width + ' which is more than 8000 permitted)';
            app.showMessage(msg, true);
            reportStatus(msg);
            return;
        }

        var params = this.getGraphicParams(pos_range[0], pos_range[1] - pos_range[0] + 1, width, true),
            formPanel = options.formPanel,
            rqCfg = params.rqCfg;
        Ext.apply(rqCfg, {
            print_title: formPanel.down('#title').checked,
            simplified: formPanel.down('#simplified').checked,
            target: options.target,
            tracks: '[key:ruler]' + params.tracks[0]
        });
        
        var processResponse = function(data, text, res) {
            if (options.formPanel.destroyed) return;

            var from_cgi = SeqView.decode(data);
            if (from_cgi.job_status) {
                if (from_cgi.job_status == 'failed' || from_cgi.job_status == 'canceled') {
                    formPanel.getEl().unmask();
                    app.showMessage(from_cgi.error_message, true)
                    reportStatus(from_cgi.error_message);
                } else {
                    if (from_cgi.progress_message && from_cgi.progress_message.length > 0) {
                        var progress_text = from_cgi.progress_message.replace(/\&quot;/gi, "\"");
                        var current_task = Ext.decode(progress_text).current_task;
                        if (current_task) {
                            app.showMessage(current_task);
                            reportStatus(current_task);
                        }
                    }
                    Ext.defer(SeqView.App.simpleAjaxRequest, 500, this, [{
                        url: params.tURLs[0] + '?job_key=' + from_cgi.job_id, context: this,
                        success: processResponse,
                        error: processResponse }]);
                }
            } else {
                var fURL = from_cgi.file_url,
                    msg = 
                formPanel.getEl().unmask();
                if (!fURL) {
                    app.showMessage('Request failed', true);
                    reportStatus('Request failed');
                    return;
                }
                // If the file_url begins with ? it contains only parameters for ncfetch, so prepend ncfetch URL
                if (fURL.charAt(0) == '?') fURL = app.m_CGIs.NetCache + fURL;

                formPanel.fileURL = fURL += '&filename=' + app.m_Config.SeqInfo.id + '[' + (rqCfg.from + 1) + '..' 
                    + (rqCfg.from + rqCfg.len) + ']' + '.' + options.target;

                reportStatus(Ext.util.Format.fileSize(from_cgi.file_size, .98)+' '+options.target + ' file created');

                if (options.nextAction == 'save') {
                    this.saveFile(formPanel.fileURL);
//                    app.showMessage('Completed');
                    if (options.target == 'pdf') formPanel.getEl().unmask();
                } else {
                    //opening
                    window.open(formPanel.fileURL + '&inline=true');
                }
            }
        }
        reportStatus('Creating File (' + options.target + ')');

        formPanel.getEl().mask('Creating File (' + options.target + ')');
        app.AjaxRequest({url: params.tURLs[0], context: this, data: rqCfg, success: processResponse, error: processResponse });
    },

    setPreferences: function() {
        var app = this.m_App;
        var optWindow = new Ext.Window({
            title: 'Preferences',
            constrain: true,
            app: app,
            layout: 'form',
            modal: true,
            buttons: [{text: 'Update',
                handler: function() {
                    var preps = {};
                    Ext.each(optWindow.items.items, function(p){preps[p.getName()] = p.getValue()});
                    Ext.apply(app.m_Config.Options, preps);
                    localStorage.setItem('NCBI/SV/Preferences', window.JSON.stringify(preps));
                    optWindow.close();
                    app.reloadAllViews();
                }},
                {text: 'Cancel', handler: function() {optWindow.close();}}]
        });

        var fields = [['color','Coloration'], [ 'label','Label Placement'], ['decor','Shapes'], ['spacing','Spacing']];
        var options = app.m_Config.Options;
        Ext.each(fields, function(p) {
            optWindow.add({
                xtype: 'combo', triggerAction: 'all', fieldLabel: p[1], mode: 'local', name: 'curr_' + p[0],
                store: options.controls[p[0]], allowBlank: false, editable: false, value: options['curr_' + p[0]]
            });
        });
        optWindow.show();
    },

    onDblClick: function(e) {
        if (!Ext.fly(e.getTarget()).hasCls('sv-dblclick')) return;
        if (this.m_Loading) return;

        var the_elem_xy = this.bodyDiv.getXY();
        var xx = e.getXY()[0] - the_elem_xy[0] - this.m_ScrollPix;

        var pix = Math.abs(xx);
        var seq_pos = this.pix2Seq(pix);

        this.zoomIn(seq_pos);
    },

    createMarkerElem: function(marker) {
        var create_params = marker.getCreateParams(false, this.m_Idx);
        return create_params.template.append(this.bodyDiv, create_params.options, true);
    },

    changeSelectedSig: function(area, ctrl_key) {
        // do nothing for "graphics" and "histograms" and for fake areas
        if (area.type === undefined || (area.type & SeqView.AreaFlags.NoSelection) !== 0)
            return;
        delete this.noParallel;
        var trk;
        if (!this.m_graphTracks.some(function(t) { trk = t; return t.id == area.parent_id; }))
            return; // return if a track not found
        var sig = area.signature + ';',
            isNew = (trk.selectedSig || '').indexOf(sig) == -1;
        if (!ctrl_key) {
            if (!isNew) isNew = trk.selectedSig != this.m_SelectedSig; //keep this one selected if there was smth else selected
            this.clearSelectedSig();
            if (isNew) {
                this.m_SelectedSig = trk.selectedSig = sig;
                delete trk.tParams;
            }
        } else {
            if (!isNew) this.clearSelectedSig(trk, sig);
            else  {
                this.m_SelectedSig = (this.m_SelectedSig || '') + sig;
                trk.selectedSig = (trk.selectedSig || '') + sig;
                delete trk.tParams;
            }
        }
        this.refresh();
      
//            this.removeSelectionsWithNoPinnedToolTips();
        this.pingClick('6-0');
        this.m_App.fireEvent('feature_clicked', this, area);
    },

    clearSelectedSig: function(track, sigList) {
        this.clearFeatMarkers();
        if (track) {
            sigList = sigList || track.selectedSig;
//            delete track.selectedSig;
        }
        if (sigList) {
            sigList.split(';').forEach(function(sig) {
                if (!sig) return;
                sig += ';'
                this.m_graphTracks.forEach(function(t) {
                    if ((t.selectedSig || '').indexOf(sig) == -1) return;
                    delete t.tParams;
                    t.selectedSig = t.selectedSig.replace(sig, '');
                    if (!t.selectedSig) delete t.selectedSig;
                    this.m_SelectedSig = this.m_SelectedSig.replace(sig, '');
                }, this);
            }, this);
        } else {
            this.m_graphTracks.forEach(function(t) {
                if (!t.selectedSig) return;
                delete t.tParams;
                delete t.selectedSig; 
            });
            this.m_SelectedSig = null;
        }
    },

    zoomRange: function() {
        var range = this.getTotalSelectedRange();
        if (range) {
            this.removeRangeSelection();
            this.startImageLoading(range[0], range[1] - range[0] + 1, {from_ui: true} );
        } else  this.gotoPositionDlg(true);
    },

    zoomTracks: function(pos, wSrc, wDst) {
        if (this.m_trackTag == 'img')
            Ext.each(this.m_graphTracks, function(t) {
                var trkTag = t.div.query('img')[0];
                if (!trkTag) return;
                trkTag.setAttribute('style', 'left:' + pos + 'px;height:' + t.imgH + 'px;width:' + wDst + 'px;position:relative;');
                t.div.setLeft(0);
            }, this);
        else
            Ext.each(this.m_graphTracks, function(t) {
                var cn = t.div.query('canvas')[0];
                if (!cn) return;
                t.div.setLeft(0);
                if (t.image && t.image.img)
                   cn.getContext('2d').drawImage(t.image.img, 0, t.imgT, wSrc, t.imgH, pos, 0, wDst, t.imgH);
            }, this);
    },

    zoomIn: function(seq_pos) {
        if (seq_pos == undefined && this.getMinLength() >= this.m_VisLenSeq) return;
        var new_len = Math.floor(this.m_VisLenSeq / 2);
        var center = seq_pos? seq_pos : this.m_VisFromSeq + new_len;
        var new_from = Math.floor(Math.max(0, center - new_len/2));
        if (new_from + new_len > this.m_App.m_SeqLength) new_len = this.m_App.m_SeqLength - new_from;
        this.zoomTracks(Math.floor((this.m_FromSeq - new_from) * this.m_bpPix * 2), this.m_Width, this.m_Width << 1);
        this.startImageLoading(new_from, new_len, {from_ui: true});
    },

    zoomOut: function(seq_pos) {
        if (seq_pos == undefined && this.m_VisLenSeq >= this.m_LenSeq) return;
        var new_len  = this.m_VisLenSeq * 2;
        var center = seq_pos? seq_pos : this.m_VisFromSeq +  this.m_VisLenSeq / 2;
        var new_from = Math.max(0, center - new_len/2);
        if (new_from + new_len > this.m_App.m_SeqLength) new_len = this.m_App.m_SeqLength - new_from;
        this.zoomTracks(Math.floor((this.m_FromSeq - new_from) * this.m_bpPix / 2), this.m_Width, this.m_Width >> 1);

        this.startImageLoading( new_from, new_len, {from_ui: true} );
    },

    zoomSeq: function(center_seq_pos) {
        var cur_center = Math.floor(this.m_VisFromSeq + this.m_VisLenSeq / 2);
        center_seq_pos = center_seq_pos || cur_center;

        var app = this.m_App,
            new_len  = this.getMinLength(),
            new_from = Math.max(0, Math.floor(center_seq_pos - new_len / 2));
        if (new_len == this.m_VisLenSeq && cur_center == center_seq_pos) return;

        if (new_from + new_len > app.m_SeqLength) new_len = app.m_SeqLength - new_from;
        var tracks_changed = false;
        if (!app.getActiveTracks().some(function(t){ return t.key == 'sequence_track' && t.shown; }))
            tracks_changed = app.showTracks({key: 'sequence_track'}, true);
        this.startImageLoading(new_from, new_len, {from_ui: true});
        if (tracks_changed) app.fireEvent('configuration_changed', app);
    },

    parseAndGotoPosition: function(text, range_only) {
        this.m_App.parseAndGotoPosition(text, {range_only: range_only, view: this});
    },

    gotoAndSearch: function(term) {
        var combo = this.m_View.down('#gotoBox');
        if (term == null || term.length == 0) {
            if (combo) {
                term = combo.getValue();
                combo.setValue('');
            }
       }
        if (!term) return;

        this.m_App.addSearchPattern(term);
        if (combo) combo.setStore(this.m_App.searchPatternData);
        this.m_App.saveSearchPatternData();

        this.m_App.gotoAndSearch(term, {view: this});
    },

    gotoPosRange: function(pos_range, center, options ){
        var new_from = pos_range[0];
        var len = this.m_VisLenSeq;
        if (pos_range.length === 2 && pos_range[1]) { // range is specified
            len =  pos_range[1] - new_from + 1;
        } else { // single position
            // according to SV-620 we want to keep the length of the current view in case if
            // it is less then 80% (?? it may change) of the sequence length
            if (len >= this.m_App.m_SeqLength * 0.8 ) {
                var new_center = new_from;
                if (new_center - 500000 < 0 && new_center + 500000 > this.m_App.m_SeqLength) {
                    if (new_center - 500 < 0 && new_center + 500 > this.m_App.m_SeqLength) {
                        var scrWidth = this.getBodyWidth();
                        new_from = new_center - scrWidth/8/2;
                        len = scrWidth/8;
                    } else {
                        new_from = new_center - 500;
                        len = 1000;
                    }
                } else {
                    new_from = new_center - 500000;
                    len = 1000000;
                }
            } else if (center) new_from -= len / 2;

            if (new_from < 0) new_from = 0;
            if (new_from + len > this.m_App.m_SeqLength) {
                new_from -= new_from + len - this.m_App.m_SeqLength;
                if (new_from < 0) {
                    new_form = 0;
                    len = this.m_App.m_SeqLength;
                }
            }
        } // end single position
        this.moveTo(new_from, len, options);
    },

    gotoPositionDlg: function(range_only) {
        Ext.MessageBox.prompt(range_only ? "Zoom to range" : "Go to position/range",
          'Enter sequence ' + (range_only ? '' : 'position or') +
              ' range<br />(possible range formats are 10k-20k, -20--10, -10k:-5, 5 to 515, -1m..1m):',
          function(btn, text) {
              if (btn == 'ok') this.parseAndGotoPosition(text, range_only);
          }, this, false);
    },

    primerBlast: function(whole, range_set) {
        //check that range is not exceeding 50Kb

        var seqinfo = this.m_App.m_ViewParams;
        if (seqinfo["acc_type"] !== "DNA") return;
        range_set = range_set || this.getSelectedRangeSet();
        if (!whole  && (!range_set || range_set.length < 1)) {
            //open message window
            Ext.MessageBox.show({
               title: 'Primer BLAST (Selection)',
               msg: 'Cannot run BLAST - select a range first',
               buttons: Ext.MessageBox.OK
           });
           return;
        } else if (!whole  && (range_set && range_set.length > 2)) {
            //open message window
            Ext.MessageBox.show({
               title: 'Primer BLAST (Selection)',
               msg: 'Cannot run Primer BLAST on more than 2 ranges',
               buttons: Ext.MessageBox.OK
           });
           return;
        } else if (range_set && range_set.length == 1){
            var range = range_set[0];
            if(range[1]-range[0]>50000) {
                Ext.MessageBox.show({
               title: 'Primer BLAST',
               msg: 'Cannot run Primer BLAST on a range longer than 50K',
               buttons: Ext.MessageBox.OK
               });
               return;
            }
        } else if (range_set && range_set.length == 2){
            var range1 = range_set[0];
            var range2 = range_set[1];
            if((range2[1]-range1[0])>50000) {
                Ext.MessageBox.show({
               title: 'Primer BLAST',
               msg: 'Cannot run Primer BLAST on ranges longer than 50K',
               buttons: Ext.MessageBox.OK
               });
               return;
            }
        }
        if (whole)
            range_set = [ [0, this.m_App.m_SeqLength-1] ];
        this.m_App.primerBlast(whole, range_set);
    },

    blastSelection: function() {
        if (!this.m_RangeSelectionSet || this.m_RangeSelectionSet.length < 1)
        {
            Ext.MessageBox.show({
               title: 'BLAST Search (Selection)',
               msg: 'Cannot run BLAST - select a range first',
               buttons: Ext.MessageBox.OK
            });
            return;
        }
        var range = this.getTotalSelectedRange();
        this.m_App.blast(this, range);
    },

    openFullView: function(mode) {
        mode = mode || 'portal';
        link = this.m_App.getLinkToThisPageURL(mode);
        if (link) {
            window.open(link + '&from=' + (this.m_VisFromSeq + 1) + '&to=' + (this.m_VisFromSeq + this.m_VisLenSeq));
            //window.open(link + (mode == 'GDV' ? ':' : '&from=') + (this.m_VisFromSeq + 1) + (mode == 'GDV' ? '-' : '&to=') + (this.m_VisFromSeq + this.m_VisLenSeq));
        }
    },

    showAssmInfo: function() {
        var sts = ['Not set in the URL (and could not be guessed)',
                'Explicitly set from the URL',
                'Automatically determined by the provided molecule ID',
                'Detection failed due to timeout',
                'Conflict between the URL parameters: molecule ID and Assembly context',
                'Automatically determined by the provided molecule ID (using cache)',
                'Molecule type known not to belong to any assembly',
                'Obtained from the molecule record'];
        var seqInfo = this.m_App.m_Config.SeqInfo;
        var msg = 'Name: <a href="' + SeqView.webNCBI + 'data-hub/genome/' + seqInfo.assm_context
            + '/" target=\"_blank\">' + seqInfo.assm_info.name
            + '</a><br>TaxID: <a href="' + SeqView.webNCBI + 'data-hub/taxonomy/' +
            + seqInfo.assm_info.tax_id + '/" target=\"_blank\">' + seqInfo.assm_info.tax_id
            + '</a><br> Type: ' + (seqInfo.assm_info.type ? 'GenBank' : 'RefSeq')
            + '<br>' + (sts[seqInfo.assm_context_status] || '');
        if (seqInfo.warnTT) msg += '<br>WARNING!<br>' + seqInfo.warnTT;
        Ext.MessageBox.show({
            icon: seqInfo.warnTT ? Ext.MessageBox.WARNING : Ext.MessageBox.INFO,
            title: 'Genomic Assembly Context',
            msg: msg,
            buttons: Ext.MessageBox.OK
        });
        return;
    },

    setGeneMode: function(mode, config) {
        var bttn = this.m_View.down('#btnGeneMode');
        if (!bttn) return;
        config = config || this.m_App.m_Config;
        if (this.m_showAllGenes != mode) config.changeGeneMode(this.m_showAllGenes = mode);
        bttn.setPressed(mode);
        bttn.setTooltip('Switch ' + (bttn.pressed ? 'OFF' : 'ON') + ' mode "show All" for Gene tracks');
        return bttn;        
    },

    hideGeneModeButton: function(hide) {
        var bttn = this.setGeneMode(hide ? false : this.m_showAllGenes || false);
        if (bttn) bttn[hide ? 'hide' : 'show']();    
    }
});  // End Of SeqView.Graphic
/*  $Id: tm_panels.js 47551 2023-08-28 15:13:40Z borodine $
*/
Ext.namespace('SeqView.TM');

Ext.define('SeqView.TM.TrackDetail', {
    extend: 'Ext.Panel',
    border: false,
    cls : 'trackDetails',
    bodyStyle: {background: 'inherit'},
    overflowY: 'auto',
    initComponent: function() {
        this.callParent(arguments);
        if (this.track) this.initTrackDetails(this.track);
    },

    initTrackDetails: function(track, title) {
        if (title) {
            this.track = track;
            this.items.each(function(item) { this.remove(item); }, this);
            this.setTitle(title);
            var help = track.help;
            if (help) {
                if (track.data_provider.indexOf('_rmt') > 0) {
                    var url = track.annots.find(function(a) { return a.search(/(https?|ftp):\/\//) == 0; });
                    if (url) 
                        help += '<br>Track source: <a href=' + url + ' target=\"_blank\">' + url.split('/').pop() + '</a><br>';
                }
                this.add({ xtype:'fieldset', padding: 5, html: help });
            }
        }

        if (!track.choice_list && !track.check_boxes) {
            this.add({xtype:'displayfield', hideLabel: true, value: "<i>No settings available</i>"});
            return;
        }
        var fPanel = this.fPanel = new Ext.FormPanel({
            border: false,
            bodyStyle: { background: 'inherit' },
            layout: 'form',
            labelWidth: 130
        });
        var self = this;
        var checkDependency = function() {
            var flds = fPanel.getForm().getFields().items;
            Ext.each(track.range_controls, function(rc) {
                self.rcForm.down('#' + rc.name).setDisabled(false);
            });
            Ext.each(flds, function(fld) {
                for (var name in fld.conditions) {
                    Ext.each(fld.conditions[name].value, function(val) {
                        if (val == fld.choiceLst.curr_value) {
                            self.rcForm.down('#' + name).setDisabled(true);
                            return false;
                        }
                    });
                }
            });
        }

        if (track.choice_list) {
            Ext.each(track.choice_list, function(lst) {
                var data = [];
                var cur_text = '';
                Ext.each(lst.values, function(v) {
                    data.push({text: (lst.prev_value == v.name ? '* ' : '') + v.display_name, value: v.name, help: v.help});
                    if (v.name == lst.curr_value) cur_text = v.display_name;
                });
                var store = Ext.create('Ext.data.Store', {
                    model: Ext.define(null, { extend: 'Ext.data.Model',
                    fields: ['text', 'value', 'help']}),
                    data: data });
                fPanel.add(Ext.create('Ext.form.field.ComboBox', {
                    allowBlank: false,
                    itemId: lst.name,
                    choiceLst: lst,
                    disabled: !track.shown,                
                    queryMode: 'local',
                    store: store,
                    editable: false,
                    fieldLabel: lst.display_name,
                    listeners: { scope: this,
                        beforeselect: function(f, r) {
                            f.choiceLst.curr_value = r.data.value; },
                        change: function(f) { if (f.conditions) checkDependency(); }
                    },
                    value: cur_text,
                    width: 250
                }));
            });
        }
        if (track.text_boxes) {
            for (var i = 0; i < track.text_boxes.length; i++) {
                var tbox = track.text_boxes[i];
                var item = new Ext.form.NumberField({
                    fieldLabel: tbox.display_name,
                    disabled: !track.shown,
                    emptyText: 'None',
                    value: tbox.value,
                    rtbox: track.text_boxes[i],
                    tooltip: tbox.help,
                    listeners: {
                        change: function(tb, nval, oval) {
                            if (nval == 0) tb.setValue('undefined');
                            tb.rtbox.value = tb.minValue == undefined ? nval : (nval || 0);
                        },
                        render: function(tb) { if (tb.value == 0) tb.setValue('undefined'); }
                        }
                    });
                if (!tbox.name == 'HistThreshold') item.setMinValue(0);
                fPanel.add(item);
            }
        }

        if (track.check_boxes) {
            var items = [];
            for(var i = 0; i < track.check_boxes.length; ++i) {
                var cbox = track.check_boxes[i];
                items.push(new Ext.form.Checkbox({
                    boxLabel: cbox.display_name,
                    checked: cbox.value,
                    disabled: !track.shown,
                    tooltip: cbox.help,
                    rcbox: cbox,
                    listeners: { change: function(cb, newVal) { cb.rcbox.value = newVal; } }
                }));
            }
            fPanel.add({
                xtype: 'checkboxgroup',
                fieldLabel: 'Other settings',
                columns: 1, //items.length < 4 ? 1 : 2,
                vertical: true,
                items: items
            });
        }
        this.add(fPanel);
        if (track.range_controls) {
//            SeqView.TM.updateRControlData(track);
            var items = [];
            for(var i = 0; i < track.range_controls.length; ++i) {
                var ctrl = track.range_controls[i];
                Ext.each(ctrl.conditions, function(c) {
                    var fld = fPanel.down('#' + c.name);
                    fld.conditions = fld.conditions || {};
                    fld.conditions[ctrl.name] = { value: c.value };
                });
                items.push({xtype: 'fieldset',
                itemId: ctrl.name,
                title: ctrl.display_name,
                checkboxToggle: true,
                checkbox: {name: 'noauto', itemId: 'noauto', //value: !ctrl.value.autoscale,
                    listeners: {change: function(cb, val) {
                        var fs = cb.ownerCt.ownerCt;
                        fs.rCtrl.value.autoscale = !val;
                        fs.down('#min').setDisabled(!val);
                        fs.down('#max').setDisabled(!val)
                    }, scope: this}
                },
                collapsible: false,
/*                listeners: {scope: this,
                    focusleave: function(rcItem) {
                        if (title && this.rcForm.getForm().isDirty()) {
                            SeqView.TM.updateRControlData(track, rcItem.rCtrl);
                        }
                    }
                },*/
                disabled: false,
                rCtrl: ctrl,
                layout: 'column',
                defaults: { columnWidth: 0.5, padding: 5, xtype: 'numberfield',
                    listeners: {
                        change: function(rc, val) { rc.ownerCt.rCtrl.value[rc.name] = val; }                    }
                },
                items: [
                    { name: 'min', itemId: 'min', fieldLabel: 'Min value', value: ctrl.value.min },
                    { name: 'max', itemId: 'max', fieldLabel: 'Max value', value: ctrl.value.max }
                ]});
            }
            this.rcForm = this.add(new Ext.FormPanel({border: false,  bodyStyle: {background: 'inherit'}, padding: 0, items: items}));
                Ext.each(track.range_controls, function(rc) {
                    var fs = this.down('#' + rc.name);
                    fs.checkboxCmp.setValue(!ctrl.value.autoscale);
                    if (Ext.isSafari12)
                        fs.on({ afterrender: function(fs) {
                            fs.titleCmp.setStyle('position', 'inherit');
                            var icb = Ext.query('div[class="x-form-cb-wrap-inner"]', fs)[0];
                            if (icb) icb.setAttribute('style', 'position:fixed;');
                        }});
                }, this);
            checkDependency();
        } else delete this.rcForm;
    }
});

SeqView.TM.updateRControlData = function(track, rctrl) {
    rctrl = rctrl || track.range_controls;
    Ext.each(rctrl, function(rc) {
        var min = rc.value.min,
            max = rc.value.max;
        if (typeof min == 'number' && typeof max == 'number' && min > max) {
            rc.value.min = max;
            rc.value.max = min;
        }
        min = parseFloat(rc.value.min);
        max = parseFloat(rc.value.max);
        rc.curr_value = (isNaN(min) ? 'inf' : min) + '|' + (isNaN(max) ? 'inf' : max) +(rc.value.autoscale ? '|auto' : '');
    });
};

SeqView.TM.collectAllTracks = function(categories, activeOnly) {
    var collected_tracks = [];
    Ext.each(categories, function(cat) {
        Ext.each(cat.subcategories, function(subcat) {
            Ext.each(subcat.tracks, function(track) {
                if (track.shown || !activeOnly ) {
                    if (cat.name && cat.name != 'other')
                        track.cat_name = cat.name;
                    else
                        delete track.cat_name;
                    if (subcat.name && subcat.name != 'other')
                        track.subcat_name = subcat.name;
                    else
                        delete track.subcat_name;
                    collected_tracks.push(track);
                }
            });
        });
    });
    collected_tracks.sort(function(t1,t2) {
        return t1.order - t2.order;
    });
    return collected_tracks;
}

SeqView.TM.collectActiveTracks = function(categories) {
    return SeqView.TM.collectAllTracks(categories, true);
}


Ext.define('SeqView.TM.TracksPanel', {
    extend: 'Ext.Panel',
    initComponent: function() {
        if (this.activeTracks) {
            this.gridPanel = new SeqView.TM.GridBase({ tracks: this.tracks, sortable: false,
                 multiSelect: true,
                 viewConfig: { plugins: { ptype: 'gridviewdragdrop', containerScroll: true },
                     listeners: { drop: function() {
                         var tracks = this.grid.tracks;
                         var order = 0;
                         this.grid.getStore().data.each(function(r) {
                              Ext.each(tracks, function(t) {
                                  if (t.cfgIdx == r.data.cfgIdx) {
                                       t.order = order++;
                                       return false;
                                  }
                              });
                         });
                     }}
                 }});
        } else 
            this.gridPanel = (this.searchTracks)
                ? SeqView.TM.createSearchTracksTab(null, this.tracks_config)
                : new SeqView.TM.GridBase({tracks: this.subcategories ? this.subcategories[0].tracks : this.tracks,
                                                        subcategories: this.subcategories});

        if (!this.searchTracks) {
            this.detailPanel = new SeqView.TM.TrackDetail({
                region: 'south',
                split: true,
                height: 225,
                minSize: 100,
                collapsible: true
            });
            var sm = this.gridPanel.getSelectionModel();
//            if (this.uploadedData) sm.setSelectionMode('MULTI');
// track details update 
           sm.on('selectionchange', this.onRowSelect, this);
        }

        var items = [];
        if (this.activeTracks) {
            var hintStr = {
                region: 'north',
                border: false,
                bodyStyle: { padding: '5px 5px 5px 5px', 'font-size': '14px' },
                html: 'Click on track to display settings. To re-order tracks, drag and drop track names.'
            };
            items.push(hintStr, this.gridPanel, this.detailPanel);
        } else if (this.searchTracks) {
            var sbar = new Ext.Panel({
                region: 'north',
                border: false,
                tbar: ['-',' Search: ', ' ', new SeqView.SearchField({store: this.gridPanel.store, width: 220, callback: this.cfgPanel.sviewApp.ping('9-6')})]
            });

            items.push(sbar, this.gridPanel);
        }
        else items.push(this.gridPanel, this.detailPanel);
        Ext.apply(this, {
            border: false,
            layout: 'border',
            items: items
        });
        this.callParent(arguments);

        this.gridPanel.on('viewready', function(g) {
            g.getSelectionModel().select(0);
        });
    },

    onRowSelect: function(sm, srow) {
        if (!srow.length) return;
        var trk = srow[0].getData();
        var title = 'Track Settings: ' + (trk.display_name || trk.name) + 
            (!trk.legend_text ? ''
            : '&nbsp&nbsp<a href="' + SeqView.getHelpURL('legends/#' + trk.legend_text)
            + '" target="_blank" style="font-weight:normal">Track legend</a>');


        this.detailPanel.initTrackDetails(trk, title);
/*        var group_button = Ext.getCmp('sv-config_group_button_id');
        var group = rdata.setting_group;
        if (group) {
            group_button.show();
            group_button.track = rdata;
            var display_name = group;
            var display_limit = 10;
            if (group.length > display_limit) {
                display_name = group.substr(0, display_limit) + "...";
            }
            group_button.setText("Configure Group '" + display_name + "'");
        } else {
            group_button.hide();
        }*/
    }
});

Ext.define('SeqView.TM.CategoryTrackPanel', {
    extend: 'Ext.Panel',
    initComponent : function(){
        Ext.apply(this, {
            title: 'Available Tracks',
            layout: 'fit',
            itemId: 'trackpanel',
            dirty: false,
            items: this.createPanel(),
            listeners: {
                show: function(p) {
                    p.updateLayout();
                }
            }
        });
        this.callParent(arguments);
    },

    reload: function(categories) {
        this.removeAll();
        this.categories = categories;
        this.dirty = false;
        var panel = this.createPanel();
        this.add(panel);
        if (this.isVisible())
            this.updateLayout();
    },

    createPanel: function() {
        var tracks_config = this.sviewApp.m_Config.TrackConfig;
        var categories = this.categories;
        var tabs = [];
        this.activeTracksTab = new SeqView.TM.TracksPanel({
            title: 'Active Tracks',
            iconCls: 'xsv-seq-logo',
            activeTracks: true,
            tracks_config: tracks_config,
            tracks: SeqView.TM.collectActiveTracks(categories),
            listeners: {
                deactivate: function() {
                    this.gridPanel.syncModifiedTracks(this.tracks_config);
                }
            }
        });
        tabs.push(this.activeTracksTab);
/*        this.searchTracksTab = new SeqView.TM.TracksPanel({
            title: 'Search Tracks',
            iconCls: 'xsv-search-button',
            searchTracks: true,
            tracks_config: tracks_config,
            cfgPanel: this,
            searchStr:'*',
            listeners: {
                deactivate: function(p) { this.gridPanel.syncModifiedTracks(); }
            }
        });
        tabs.push(this.searchTracksTab);*/
        function updateRButton(sm) {
            Ext.getCmp('sv-delete_button_id').setText('Remove track(s): ' + sm.getCount() + '/' + sm.store.count())
            .setTooltip('Use (Ctrl or Shift)+click to select group of tracks for removal<br>Selected: ' + sm.getCount() + '/' + sm.store.count());
        }
        for (var i = 0; i < categories.length; i++) {
            var cat = categories[i];
            var noUserData = false; // cat.name.search(/Uploaded|Remote/) == 0;
            Ext.each(cat.subcategories, function() {
                Ext.each(this.tracks, function(){ return noUserData = this.id.search(/U|R/) != 0; });
                return noUserData;
            });
            var catTab = new SeqView.TM.TracksPanel({
                title: cat.display_name || cat.name,
                layout: 'fit',
                category: cat,
                subcategories: cat.subcategories,
                tracks_config: tracks_config,
                noUserData: noUserData,
                listeners: {
                    deactivate: function(p) {
                        Ext.getCmp('sv-delete_button_id').disable()
                        .setText('Remove track(s)')
                        .setTooltip('Remove track(s) from NCBI data storage'); 
                        this.gridPanel.syncModifiedTracks(this.tracks_config);
                    }
                }
            });
            if (!noUserData) {
                catTab.addListener('activate', function(p) {
                        Ext.getCmp('sv-delete_button_id').enable();
                        updateRButton(p.gridPanel.getSelectionModel());
                    });
                catTab.gridPanel.addListener('selectionchange', updateRButton);
            }
            tabs.push(catTab);
        }
        var panel = new Ext.tab.Panel({
            activeTab: 0,
            tabPosition: 'left',
            tabRotation: 0,
            plain:true,
            tabWidth: 158,
            items: tabs
        });
        return panel;
    },
    
    getActiveGridPanel: function() {
        return this.items.items[0].getActiveTab().gridPanel;  
    },

    setDirty: function(flag) {
        this.dirty = (flag === true);
    }
});

Ext.define('SeqView.TM.UploadPanel', {
    extend: 'Ext.Panel',
    initComponent : function(){
        Ext.apply(this, {
            title: 'Custom Data',
            layout: 'border',
            itemId: 'uploadpanel',
            infoMsg: 'No informations or details',
            items:[{
                region: 'west', layout:'fit', width: 168, title:'Data Source',
                items: [{
                    xtype: 'treepanel', itemId: 'treepanel', rootVisible: false, enableDD: false, lines: false,  autoScroll: false,
                    listeners: { itemClick: function(node, rec) {
                        this.uploadButt.disable();
                        this.updateDiffPart(rec.data.type, rec.data.text);
                    }, scope: this},
                      root: {
                        children: [
                            {text:'BLAST Results', type:'rid', iconCls:'xsv-ext_data', leaf:true},
                            {text:'Data File', type:'file', iconCls:'xsv-ext_data', leaf:true},
//                            {text:'BAM Files', type:'BAMfiles', iconCls:'xsv-ext_data', leaf:true},
                            {text:'Remote Files', type:'BAMremote', iconCls:'xsv-ext_data', leaf:true},
                            {text:'Alignment MUSCLE/FASTA', type:'FASTAfile', iconCls:'xsv-ext_data', leaf:true},
                            {text:'URL', type:'url', iconCls:'xsv-ext_data', leaf:true},
                            {text:'Text', type:'text', iconCls:'xsv-ext_data', leaf:true}
                    ]}
                }]
            },{
                frame:true,
                xtype:'form',
                region: 'center',
                hidden: false,
                itemId:'loaddatapanel',
                fileUpload: true,
                items:[                        
                    {xtype:'panel', layout: 'form', border: false, bodyStyle: {background: 'inherit'},
                    items:[
                        {xtype:'textfield', fieldLabel: 'Track Name', name:'track_name', emptyText: 'Optional parameter'},
                        {xtype:'hidden', name:'assm_acc', value: this.cfgPanel.sviewApp.m_AssmContext || ""}]},
                    {xtype:'displayfield', id: 'sv-uplmsg' + this.cfgPanel.sviewApp.m_Idx, value: '', width:'95%',
                        style: {color:'grey', 'white-space': 'nowrap'}},
                    {xtype:'displayfield', id: 'sv-uplerr' + this.cfgPanel.sviewApp.m_Idx, value: '',
                        style: {color:'red'}},
                    {xtype: 'button', hidden: true, scope: this, handler: this.showErrorDetails,
                        id: 'sv-err_details_button_id'}
                ]
            }],
            buttons: [{
                text:'Upload', scope: this, itemId: 'uploadButton',
                disabled: true,
                handler: function(b, e) {
                    e.stopEvent();
                    var form = this.down('#loaddatapanel').getForm();
                    if (!form.isValid()) return;
                    if (form.isDirty() || this.dzfiles) {
                        this.submitData(form);
                        var data = {'sv-area': '9-5-' + this.type, jsevent: 'upload', name: form.getValues().track_name};
                        if (!data.name) delete data.name;
                        this.cfgPanel.sviewApp.ping(data);
                    }
                    else {
                        this.cfgPanel.sviewApp.showMessage('Inputs are empty or invalid!', true);
                    }
                }
            }]
        });
        this.callParent(arguments);
        this.down('#treepanel').getSelectionModel().select(0);
        this.uploadButt = this.down('#uploadButton');
        this.updateDiffPart('rid', 'BLAST Results');
    },
    showErrorDetails: function() {
        if (this.infoMsg == '') return;
        var mbw = new Ext.Window({
            title: 'Error details',
            constrain: true,
            modal: true, layout:'fit',
            width: 400, height: 300});
        mbw.add(new Ext.FormPanel({ labelWidth: 1, autoScroll: true,
            items:[
               {xtype: 'displayfield', value: this.infoMsg, textalign: 'left'}],
        	   buttons:[{text: 'OK', handler: function() {mbw.close();}}]})
        );
        mbw.show();
    },
    Sec2Time: function(sec) {
        var hours   = Math.floor(sec / 3600);
        var minutes = Math.floor((sec - (hours * 3600)) / 60);
        var seconds = sec - (hours * 3600) - (minutes * 60);
        var time = '';
        if (hours > 0) 
            time = hours + ' hour' + ((hours > 1) ? 's ' : ' ');
        if (minutes > 0)
            time += minutes + ' minute' + ((minutes > 1) ? 's ' : ' ');
        if (seconds > 0)
            time += seconds + ' second' + ((seconds > 1) ? 's' : '');
        return time;
    },
    updateMessage: function() {
        this.currTime++;
        this.totalTime++;

        var msg = '<br>';
        msg += 'Total time: ' + this.Sec2Time(this.totalTime) + '<br>';
        for (var i = 0; i < this.tasks.length; i++)
            msg += this.tasks[i].task + ": " + this.tasks[i].time + " seconds<br>";
        msg += this.currTask + ": " + this.Sec2Time(this.currTime);
        if (this.percentage) {
            var time = this.Sec2Time(Math.round(this.currTime * (100/this.percentage - 1)));
            if (time.length > 0) 
                msg += " (" + time + " remaining)";
        }
        msg += '<br>';
        if (this.percentage) msg += 'Percentage: ' + this.percentage + '%<br>';
        this.cfgPanel.sviewApp.showMessage(msg);
        Ext.defer(this.updateMessageWrap, 1000, this, []);
    },

    updateMessageWrap: function() {
        if (this.uploaderUUD) this.updateMessage("");
    },

    reset: function(panel) {
        if (this.dzfiles) {
            for (var f in this.dzfiles) {
                this.dzfiles[f].dzfname.update('');
                this.dzfiles[f].inputFld.value = '';
            }
            delete this.dzfiles;
        } else 
            panel.getForm().setValues({rid:'', data: '', dataURL: ''});
    },

    cleanupUpload: function(msg) {
        delete this.uploaderUUD;

        var button = Ext.getCmp('sv-configure_button_id')
        if (button) button.enable();
        if (msg) this.cfgPanel.sviewApp.showMessage(msg, true);
    },
    
    submitData: function(form, callback, scope) {
        this.currTime = this.totalTime = this.percentage = 0;
        this.currTask = 'uploading file';
        this.tasks = [];

        Ext.getCmp('sv-err_details_button_id').hide();
        Ext.getCmp('sv-configure_button_id').disable();
        this.uploadButt.disable();
        this.cfgPanel.sviewApp.showMessage('');
        this.cfgPanel.sviewApp.showMessage('', true);

        var config = {assm_acc: this.cfgPanel.sviewApp.m_AssmContext};
        var fval = form.getValues();
        for (var v in fval) {
            var val = fval[v];
            if (typeof val == 'string' && val.length == 0) continue;
            config[v] = v == 'dataURL' ? val.trim() : val;
        }
        var find_comp = config.find_comp || false;
        delete config.find_comp;
        if (config.rid) config.blast = {rid: config.rid, link_related_hits: find_comp};
        else config.check_cs = true;
        var uPanel = this;
        if (callback) {
            this.exitCallback = {callback: callback, scope: scope};
        }
        uPanel.consError = console.error;
        console.error = function() {}
        var finalize = function(msg) {
            console.error = uPanel.consError || console.error;
            delete uPanel.consError;
            uPanel.cleanupUpload(msg);
            if (!msg) {
                var msg = 'Data uploaded/registered';
                if (uPanel.cfgPanel.sviewApp.m_AssmContext) {
                    msg += " on assembly " + uPanel.cfgPanel.sviewApp.m_AssmContext;
                }
                uPanel.cfgPanel.sviewApp.showMessage(msg);
                var dTrk = Ext.each(uPanel.infoUpTracks, function() { return (this.track_type == 'non-displayable_track') });
                if (typeof dTrk == 'undefined') uPanel.cfgPanel.sviewApp.showMessage('There is no data that can be displayed on the sequence', true);
            }
            uPanel.reset(uPanel.down('#loaddatapanel'));
            if (uPanel.exitCallback) {
                var tmpobj = uPanel.exitCallback;
                delete uPanel.exitCallback;
                tmpobj.callback.call(tmpobj.scope);
            }
        }

        if (this.dzfiles) Ext.apply(config, this.dzfiles);
        try {
            this.uploaderUUD = SeqView.App.getFileUploader(config);
            var promise = this.uploaderUUD.getPromise();
            promise.fail(function(){
                finalize('Failed to process data: ' + this.getErrors());
            });
            promise.done(function(tlist, dkey) {
                uPanel.infoUpTracks = this.getTracks();
                if (uPanel.infoUpTracks) {
                    Ext.getCmp('sv-load_default_button_id').disable();
                    uPanel.cfgPanel.sviewApp.addUUDTracks(uPanel.infoUpTracks);
                    uPanel.ownerCt.down('#trackpanel').setDirty(true);
                }
                var errMsg = this.getErrors();
                if (errMsg.length) {
                    uPanel.infoMsg = '';
                    Ext.each(errMsg, function(msg, idx){
                        this.infoMsg +='# ' + (idx + 1) + '. ' + msg + '<br>';
                    }, uPanel);
                    
                    var bttn = Ext.getCmp('sv-err_details_button_id');
                    var bttxt ='Data parsing error details (' + errMsg.length + ')';
                    bttn.setText(bttxt);
                    bttn.show();
                }
                finalize();
            });
            promise.progress(function(progress) {
                if (progress.percentage) uPanel.percentage = progress.percentage;
                var task = progress.current_task;
                if (task == "" || task == 'pending' || task == uPanel.currTask) return;

                uPanel.tasks.push({task: uPanel.currTask, time: uPanel.currTime});
                uPanel.currTask = task;
                uPanel.currTime = 0;
            });
            this.uploaderUUD.upload();
        } catch(e) { finalize('Unable to process data: ' + e.message); } 
    
        this.updateMessage();
    },

    updateDiffPart: function(type, title) {
        var uPanel = this;
        this.type = type;
        var fileUploadDisabled = Ext.isIE && !(Ext.isIE11 || Ext.isIE10);
        var manageButt = function(str) {
            uPanel.uploadButt[str.length ? 'enable' : 'disable']();
            return true;
        };
        var processFName = function(ffld, file) {
            file.inputFld = document.getElementById(ffld.getInputId());
            file.dzfname = uPanel.down('#dz_' + ffld.itemId);
            uPanel.dzfiles = uPanel.dzfiles || {};
            uPanel.dzfiles[ffld.itemId] = file;
            file.dzfname.update(file.inputFld.value = file.name);
        }
        
        var updateDropZone = function(self) {
            var files = self.fileInputEl.dom.files;
            if (files.length == 0) {
                uPanel.down('#dz_' + self.itemId).update('');
                uPanel.uploadButt.disable();
            } else {
                processFName(self, files[0]);
                uPanel.uploadButt.enable();
            }
        };
        var helpTxt = fileUploadDisabled ? 'File upload is unavailable in Internet Explorer versions 9 and earlier' : '';
        if (SeqView.jsonType == 'JSONP') {
            fileUploadDisabled = true;
            helpTxt = 'Local files uploading is currently unavailable for X-domain/IE configuration';
        }
        var diff_parts = {
            rid: [
                {xtype:'component', padding: '5px', html: 'Please enter NCBI BLAST request ticket (RID) then press Upload to add new alignment track.'},
                {xtype:'textfield', hideLabel:true, emptyText:'Please enter Blast RID', name: 'rid', width:'99%',
                     validator: manageButt},
                {xtype:'checkbox', name: 'find_comp', height:15, boxLabel:'Link related hits together',  checked:true, style: {fontSize: '100%'}},
                {xtype:'component', padding: '5px', html: 'BLAST returns separate alignments for each query, and these separate alignments can further be ordered into sets offering consistent non-overlapping query and subject coverage. The sequence viewer offers the ability to evaluate the original BLAST hits on-the-fly and link together alignments that meet a strict definition of non-overlapping query and subject coverage.'}
            ],
            text: [
                {xtype:'component', padding: '5px', html: 'Please paste you text annotations then press Upload to add new track.'},
                {xtype:'textarea', hideLabel:true, emptyText:'Please paste your text here', name: 'data', width:'99%', height:135,
                     validator: manageButt}
            ],
            url: [
                {xtype:'component', padding: '5px', html: 'component', padding: '5px', html: 'Please specify the WEB URL to download the input file then press Upload to add new track(s).'},
                {xtype:'textfield', hideLabel:true, emptyText:'Please enter URL', name: 'dataURL', width:'99%', validator: manageButt}
            ],
            BAMremote: [
                {xtype:'component', padding: '5px', html: 'Specify the URL of a remote file, then press Register to add the new track. Currently, BAM, BigBed, BigWig, VCF, and indexed VCF data types are supported. BAM files must have a .bai index file at the same location.'},
                {xtype:'hidden', name:'remote', value: true},
                {xtype:'hidden', name:'hint_acc', value: this.cfgPanel.sviewApp.GI},
                {xtype:'textfield', hideLabel:true, emptyText:'Please enter URL', name: 'dataURL', width:'99%', validator: manageButt}
            ],
            BAMfiles: (helpTxt)
                ? [{xtype:'component', padding: '5px', html:  helpTxt}]
                : [
                   {xtype:'hidden', name:'dst', value: 'filetrack'},
                   {xtype:'hidden', name:'file_format', value: 'bam'},
                   {xtype:'component', padding: '5px', html: 'Please specify or drop input files then press Upload to add new track(s).'},
                   {xtype:'panel', layout: 'form',  border: false, bodyStyle: {background: 'inherit'},
                        items: [
                            {xtype: 'filefield',  itemId: 'file', fieldLabel: 'File to upload', listeners: {change: updateDropZone}},
                            {xtype: 'filefield',  itemId: 'auxiliary_file', fieldLabel: 'Index file', listeners: {change: updateDropZone}}
                        ]
                   },
                   {xtype: 'fieldset', title: 'Drag and drop file here', id: 'sv_dropzone', name: 'bam', height: 60, labelWidth: 1, width: '100%',
                        items: [{xtype: 'displayfield', itemId: 'dz_file', value: ''},
                                {xtype: 'displayfield', itemId: 'dz_auxiliary_file', value: ''}]}
                ],
            FASTAfile: (helpTxt)
                ? [{xtype:'component', padding: '5px', html: helpTxt, style: {color: 'red'}}]
                : [{xtype:'component', padding: '5px', html: 'Please specify or drop input MUSCLE/FASTA file then press Upload to add new track(s).'},
                   {xtype:'hidden', name:'file_format', value: 'align'},
                   {xtype:'panel', layout: 'form',  border: false, bodyStyle: {background: 'inherit'},
                        items: [
                            {xtype: 'filefield',  itemId: 'file', fieldLabel: 'File to upload', listeners: {change: updateDropZone}}
                        ]
                   },
                   {xtype: 'fieldset', title: 'Drag and drop file here', id: 'sv_dropzone', name: 'align', height: 60, labelWidth: 1, width: '100%',
                        items: [{xtype: 'displayfield', itemId: 'dz_file', value: ''}]}
                ],

            file: (helpTxt)
                ? [{xtype:'component', padding: '5px', html: helpTxt, style: {color: 'red'}}]
                : [{xtype:'component', padding: '5px', html: 'Please specify or drop an input file then press Upload to add new track(s).'},
                   {xtype:'panel', layout: 'form',  border: false, bodyStyle: {background: 'inherit'},
                        items: [
                            {xtype: 'filefield',  itemId: 'file', fieldLabel: 'File to upload', listeners: {change: updateDropZone}},
                            {xtype: 'combo', triggerAction:'all', width:120, name:'file_format', fieldLabel:'File format', mode:'local',
                                store:['auto detect', 'asn binary', 'asn text', 'bed', 'bed 15', 'fasta', 'gff2', 'gff3', 'gtf/gff', 'gvf', 'vcf', 'wiggle', 'xml'],
                                allowBlank: false, editable: false, value: 'auto detect'}
                        ]
                   },
                   {xtype: 'fieldset', title: 'Drag and drop file here', id: 'sv_dropzone', height: 60, labelWidth: 1, width: '100%',
                        items: [{xtype: 'displayfield', itemId: 'dz_file', value: ''}]}
                ]
            };
        var panel = this.down('#loaddatapanel');
        var dp = panel.down('#diff_part');
        if (dp) panel.remove(dp);
        panel.insert(0, new Ext.form.Panel({items: diff_parts[type], itemId: 'diff_part', border: false, bodyStyle: {background: 'inherit'}}));
        this.uploadButt.setText(type.search(/remote/) < 0 ? 'Upload' : 'Register');
        panel.setTitle(title);
        panel.updateLayout();

        this.cfgPanel.sviewApp.showMessage(type.search(/file|url/) >= 0 ? 'The current per-file limit for uploads is 4GB' : '');

        if (!fileUploadDisabled && type.indexOf('file') >= 0) {   
            var dzone = document.getElementById('sv_dropzone');
            dzone.ondragover = function(){ return false; };
            dzone.ondragleave = function(){ return false; };
            var dtfile = uPanel.down('#file'),
                idxfile = uPanel.down('#auxiliary_file');
            dzone.ondrop = function(e) {
                e.preventDefault();
                var files = e.dataTransfer.files;
                if (type == 'BAMfiles') {
                    Ext.each(files, function(f, idx) {
                        processFName(/.(BAM)$/i.test(f.name) ? dtfile : idxfile, f);
                    }, uPanel, true);
                    if (!uPanel.dzfiles || !uPanel.dzfiles.file || !uPanel.dzfiles.auxiliary_file) return;
                }
                else processFName(dtfile, files[0]);
                uPanel.uploadButt.enable();
                uPanel.cfgPanel.sviewApp.pingClick('9-5-DnD');
                return false;
            };
        }
    }
});

Ext.define('SeqView.tmsTrackSearch', {
    extend: 'Ext.data.Model',
    fields: [ 'group', 'name', 'description', 'shownBool'],
    proxy: {
        type: 'ajax',
        url: 'tms/search',
        reader: {
            rootProperty: 'tracks',
            totalProperty: 'count'
        },
        // sends single sort as multi parameter
        simpleSortMode: true
    }
});


Ext.define('SeqView.TM.gridSearchPaging', {
    extend: 'Ext.grid.Panel',
    xtype: 'paging-grid',
    multiSelect: true,
    itemId: 'pgGrid',
    height: '100%',
    flex: 1,
    loadMask: true,
    
    initComponent: function(){
        var trackFilter = function(rec) {
            if (!this.baseParams.regexp) return true;
            try {
                if (rec.data.help.search(this.baseParams.regexp) >= 0) return true;
                if (rec.data.display_name.search(this.baseParams.regexp) >= 0) return true;
            } catch(err) { console.log(err.message); }        
            return false;
        };
        var isTMS = this.tracks == false,
            store = new Ext.data.Store({
                data: this.tracks,
                model: isTMS ? 'SeqView.tmsTrackSearch' : 'SeqView.TrackList',
                pageSize: 100,
                baseParams: {query: this.srchStr || '', limit: 20, filterBy: trackFilter, regexp: this.srchStr ? globStringToRegex(this.srchStr) : ''},
                remoteSort: true,
                sorters: [{
                    property: 'Category',
                    direction: 'ASC'
                }]
            });
        if (!isTMS) this.tracks.forEach(function(t) { delete t.subcat; });

        function setTip(value, metaData) {
            metaData.tdAttr = 'data-qtip="' + value + '"';
            return value + (Ext.isSafari ? '<div></div>' : '');
        }

        Ext.apply(this, {
            store: store,
            columns:[
                {text: "Category", width: 70, sortable: !isTMS, dataIndex: (isTMS ? 'group' : 'subcat'), editable: false, menuDisabled: true},
                {text: 'Active', width: 50, xtype: 'checkcolumn', sortable: !isTMS, dataIndex: 'shown' + (isTMS ? 'Bool' : ''), menuDisabled: true,
                    listeners: {checkchange: function(col, num, tick, trk) {
                        var sel = this.getSelection();
                        if (sel.find(function(t){ return t == trk; }))
                            sel.forEach(function(t){ t.set(col.dataIndex, tick); this.updateTrkList(t, tick);  }, this);
                        else this.updateTrkList(trk, tick);
                    }, beforecheckchange: function(col, num, tick, trk, opts) {
                       var ret = true;
                       return ret;
                    }, scope: this},
                    editor: { xtype: 'checkbox', cls: 'x-grid-checkheader-editor'}},
                {text: "Track name", sortable: true, dataIndex: 'display_name', flex: true, editable: false, menuDisabled: true, renderer : setTip},
                {text: "Description", dataIndex: 'help', flex: true, editable: false, menuDisabled: true, renderer: setTip}

            ]
        });
        if (isTMS) this.bbar = {
                xtype: 'pagingtoolbar',
                listeners: {change: function(ptb, pgpar) {
                    var trx = this.store.data.items;
                    this.checkedTracks.forEach(function(t, idx, arr) {
                        if (!t.modified) return;
                        var trk = trx.find(function(tr) { return t.id == tr.id; });
                        if (!trk) return;
                        trk.set('shownBool', Boolean(t.data.shown));
                        trk.data.shown = t.data.shown;
                        arr[idx] = trk;
                    }
                    );
                }, scope: this},
                displayInfo: true,
                displayMsg: 'Displaying tracks {0} - {1} of {2}',
                emptyMsg: "No tracks to display"
            }
        else this.bbar = ['->', 'Tracks available: ', {itemId: 'foundTracks', xtype: 'tbtext', text: ''}];
        this.callParent(arguments);
        if (this.srchStr) store.filterBy(trackFilter, store);
    },
    
    afterRender: function(){
        this.callParent(arguments);
        if (this.tracks) {
            this.addListener('filterchange', this.updateBBarInfo, this);
            this.getDockedItems('toolbar')[0].setHeight(30);
            this.updateBBarInfo();
        }
        this.getStore().loadPage(1);
    },

    setExtraParams: function(params) {
        this.store.getModel().getProxy().setExtraParams(params);
        if (params.track_class_scs[0] != 'Phenotype*') this.columns[2].renderer =
            function setTipWithLink(value, metaData, rec) {
                metaData.tdAttr = 'data-qtip="' + value + '"';
                var smp = rec.data['GEO Sample'] || value.split(':', 1)[0];
                    url = rec.data.url || (SeqView.webNCBI + 'sra\/' + rec.data.SRARun);
                return value.replace(smp, '<a target="_blank" href="' + url + '">' + smp + '</a>') + (Ext.isSafari ? '<div></div>' : '');
            }

    },    

    setTerm: function(term) {
        this.store.getModel().getProxy().setExtraParam('t', term);
    },

    updateBBarInfo: function() {
        this.down('#foundTracks').setText(this.store.getData().length + '/' + this.tracks.length);
    },
    
    updateTrkList: function(trk, tick) {
        var idx = this.checkedTracks.findIndex(function(t) { return t.id == trk.id; });
        if (idx >= 0) this.checkedTracks[idx] = trk;
        else if (tick) this.checkedTracks.push(trk);
        trk.data.shown = tick;
    }
});


Ext.define('SeqView.TM.SearchPanel', {
    extend: 'Ext.Panel',
    title: 'Find Tracks',
    itemId: 'searchPanel',
    iconCls: 'xsv-search-button',
    layout: 'box',
    initComponent: function() {
        var app = this.cfgPanel.sviewApp,
            srchPan = this,
            seqInfo = app.m_Config.SeqInfo,
            srchFld = new SeqView.SearchField({width: 160, callback: function() {
                if (this.ownerCt.down('#availTrx').getValue()) app.pingClick('9-6-avail');
            }});
        this.checkedTracks = SeqView.TM.collectActiveTracks(this.categories);
        if (!SeqView.extCreateRequest) {
            SeqView.extCreateRequest = Ext.Ajax.createRequest;
            Ext.Ajax.createRequest = function(opts, rqOpts) {
                if (rqOpts.url.indexOf('tms/search') == 0) {
                    var ret = SeqView.extCreateRequest.call(this, opts, rqOpts);
                    var params = {},
                        optSet = rqOpts.url.split('&');
                    optSet.shift();
                    optSet.forEach(function(opt) {
                        opt = opt.split('=');
                        switch (opt[0]) {
                         case 'sort':
                          if (opt[1]=='help')
                          {
                            params[opt[0]]="description_s"; break;
                            }
                          if (opt[1]=='display_name')
                            params[opt[0]]="track_name_s"; break;


                            case 'dir': case 'page': case 'resolve': case 't': case 'assm': params[opt[0]]=opt[1]; break;
                            case 'limit': params.count = opt[1]; break;
                            case 'track_class_scs': params[opt[0]]=[opt[1]]; break;
                        }
                    });
                    rqOpts.url = app.m_CGIs.FeatSearch + '?version';
                    rqOpts.url.replace('search', 'watchdog');
                    
                    var callback = opts.callback;
                    opts.callback = function(o, s, r) {
                        var scope = this,
                            chTrx = srchPan.checkedTracks;//.filter(function(t) { return t.dirty; });
                        SeqView.tmsSearch(params, function(data) {
                            data.tracks.forEach(function(t) {
                                var idx = chTrx.findIndex(function(tt) { return tt.id == t.id; });
                                if (idx < 0) {
                                    t.shownBool = false;
                                    return;
                                }
                                if (!chTrx[idx].internalId) { t.shown = chTrx[idx].shown; chTrx[idx].modified = true; chTrx[idx].data = {shown: t.shown}; }
                                var trk = chTrx[idx];
                                t.shown = t.shownBool = typeof trk.modified == 'object' && trk.modified.hasOwnProperty('shownBool')
                                    ? trk.modified.shownBool : Boolean(trk.data.shown);
                            });
                            r.responseJson = data;
                            callback.call(scope, o, s, r);
                        });
                    };
                    return ret;
                }
                return SeqView.extCreateRequest.call(this, opts, rqOpts);   
            }
        }

        var mp = new Ext.Panel({
            border: false,
            cls: 'SeqViewerApp',
            height: '100%', width: '100%',
            layout: 'hbox'
        });
        mp.add(
            {height: '100%', width: 185,
             items  :[{
                xtype: 'fieldset',
                height: '100%', width: '100%',
                border: 'none',
                title: 'Search',
                items: [ srchFld,
                    {xtype: 'radiogroup', columns: 1, vertical: true, itemId: 'rbtns',
                     defaults: {name: 'searchRB', listeners: {scope: this,
                        change: function(rf, val) {
                            if (!val) return;
                            mp.remove(mp.down('#pgGrid'), false);
                            rf.gridPanel = rf.gridPanel || new SeqView.TM.gridSearchPaging({ title: rf.boxLabel, checkedTracks: this.checkedTracks,
                                 tracks: rf.inputValue == 'Avail' ? app.m_Config.TrackConfig : false });
                            switch (rf.inputValue) {
                                case 'Avail':
                                    delete srchFld.userClickTrigger;
                                    break;
                                default: 
                                    rf.gridPanel.setExtraParams({assm: seqInfo.assm_context, resolve: 1, track_class_scs: [rf.inputValue]});
                                    srchFld.userClickTrigger = function(term) {
                                        rf.gridPanel.setTerm(term);
                                        rf.gridPanel.store.removeAll();
                                        rf.gridPanel.store.loadPage(1);
                                    }
                            }
                            mp.add(rf.gridPanel);
                        }}
                    },
                    items: [
                        {boxLabel: 'Available Tracks', inputValue: 'Avail', itemId: 'availTrx'},
                        {boxLabel: 'GEO samples', inputValue: 'GEO*'},
                        {boxLabel: 'SRA alignments', inputValue: 'SRA*', disabled: seqInfo.assm_context == ''},
                        {boxLabel: 'dbGaP aligned analyses', inputValue: 'Phenotype*'}
                    ]
                }]
            }]
            });
        this.callParent(arguments);
        this.add(mp);
        this.on('beforedestroy', function() {
            Ext.Ajax.createRequest = SeqView.extCreateRequest || Ext.Ajax.createRequest; delete SeqView.extCreateRequest;
        });

        mp.down('#rbtns').setValue({searchRB:'Avail'});
        this.availTracksStore = srchFld.store = mp.down('#availTrx').gridPanel.getStore();
    }
});


/*  $Id: tm_dialog.js 46965 2022-03-16 19:23:10Z evgeniev $
*/

function encodeKeyVal(s, key, val, f_escape) {
    if (val && val.length > 0) {
        s += (s.length > 0 ? ',' : '') + key + ':' + (f_escape ?  NCBIGBUtils.escapeTrackName(val) : val);
    }
    return s;
}

SeqView.TM.isDownloadable = function(track) {
    var dKeys = 'gene_model_track,aggregate_feature_track,SNP_track';
    if (track.data_provider != 'seqgraphic_rmt')
        return dKeys.indexOf(track.key) >= 0;
    else
        return ((track.dbname=="vcfTabix") || (track.dbname=="vcf"));

}

SeqView.TM.processTracksInfo = function(cgi_tracks) {
    var ret = [];
    var search_category = function(categories, cat_name) {
        var scat = null;
        Ext.each(categories, function(cat) {
            if (cat.name == cat_name) {
                scat = cat;
                return false;
            }
        });
        return scat;
    }
    
    Ext.each(cgi_tracks, function(ctrack) {
        var cat = null;
        ctrack.name = ctrack.name || ('track_' + ctrack.display_name || ctrack.cfgIdx);
        ctrack.display_name = ctrack.display_name || ('Track ' + ctrack.cfgIdx);
        if (!ctrack.order && ctrack.order !== 0) ctrack.order = cgi_tracks.length;
        var otherCat = {name: 'other', display_name: 'Other', help: 'Other tracks', order: 999999};
        if(ctrack.category) {
            cat = search_category(ret, ctrack.category.name)
            if (!cat) {
                cat = Ext.clone(ctrack.category);
                cat.subcategories = [];
                ret.push(cat);
            }
        } else {
            cat = search_category(ret, 'other')
            if (!cat) {
                cat = otherCat;
                cat.subcategories = [];
                ret.push(cat);
            }
        }
        var subcat = null;
        if (ctrack.subcategory) {
            subcat = search_category(cat.subcategories, ctrack.subcategory.name)
            if (!subcat) {
                subcat = Ext.clone(ctrack.subcategory);
                subcat.tracks = [];
                cat.subcategories.push(subcat);
            }
        } else {
            subcat = search_category(cat.subcategories, 'other')
            if (!subcat) {
                subcat = otherCat;
                subcat.tracks = [];
                cat.subcategories.push(subcat);
            }
        }
        subcat.tracks.push(ctrack);
    });
    ret.sort(function(a, b) { return a.order - b.order; });
    return ret;
};

////////////////////////////////////////////////////////////////////////////////////////////////
//
SeqView.TM.trackToString = function(track, for_group, tmsSearch) {
     if (track.tParams && track.tParams.charAt(0) !== '[') return track.tParams;

    var str = '';

    str = encodeKeyVal(str, 'key', track.key);
    if (for_group) {
        str = encodeKeyVal(str, 'group_name', track.setting_group);
    } else {
        Ext.each(['name', 'display_name', 'id', 'data_key', 'subkey', 'dbname','big_data_index'], 
            function(p) { str = encodeKeyVal(str, p , track[p], true);
        });
        str = encodeKeyVal(str, 'category', track.cat_name);
        str = encodeKeyVal(str, 'subcategory', track.subcat_name);
        str = encodeKeyVal(str, 'uid', track.uId, true);
        str = encodeKeyVal(str, 'setting_group', track.setting_group);
        Ext.each(track.subTracks, function(subtrack, idx) {
            str += idx == 0 ? ',subtracks:' : '|';
            str += NCBIGBUtils.escapeTrackName(subtrack);
        });
        Ext.each(track.annots, function(annot, idx) {
            str += idx == 0 ? ',annots:' : '|';
            str += NCBIGBUtils.escapeTrackName(annot);
        });
        Ext.each(track.comments, function(comment, idx) {
            str += idx == 0 ? ',comments:' : '|';
            str += NCBIGBUtils.escapeTrackName(comment.label) + '|' + comment.pos_str;
        });
        Ext.each(track.highlights, function(hl, idx) {
            str += idx == 0 ? ',highlights:' : '|';
            str += NCBIGBUtils.escapeTrackName(hl);
        });
        if (tmsSearch && track.tmsSearch) str += ',help:' + NCBIGBUtils.escapeTrackName(track.help);
        if (track.filter && track.filter.length > 0) {
            // Do not escape the bar '|'
            str += ',filter:' + NCBIGBUtils.escapeTrackName(track.filter, /([\]\[\\,:=&;"#%])/);
        }
    }
    str += SeqView.TM.getTrackDisplayOptions(track);
    Ext.each(['highlights_color', 'is_private', 'batch', 'stored_scale', 'extra_info', 'rmt_mapped_id'], 
        function(p) { if (track[p]) str += ',' + p + ':' + track[p]; });
    return str;
}

SeqView.TM.getTrackDisplayOptions = function(track) {
    var str = '';
    SeqView.TM.updateRControlData(track);
    Ext.each(['choice_list', 'check_boxes', 'hidden_settings', 'text_boxes', 'range_controls'], function() {
        Ext.each(track[this], function() {
            str += ',' + this.name + ':' + (this.hasOwnProperty('curr_value') ? this.curr_value : this.value);
        });
    });
    return str;
}

SeqView.TM.tracksArrayToString = function(tracks, shown, order, extraOpt) {
    extraOpt = extraOpt || {};
    var str = '',
        extraStr = extraOpt.extraStr || '';
    Ext.each(tracks, function(track) {
        str += '[' + SeqView.TM.trackToString(track, null, extraOpt.tmsSearch);
        if (shown) str += ',shown:' + (track.shown ? 'true' : 'false');
        if (order) str += ',order:' + track.order;
        str += extraStr + ']';
    });
    return str;
}

// returns an array of strings for parallel rendering
SeqView.TM.tracksToArrayOfStrings = function(tracks, extraOpt){
    if (!tracks.length) return {tURLs:[]};
    var trx = [],
        tURLs = [],
        affTokens = [];
    for (var i = 0; i < tracks.length; i++) {
        if (!tracks[i].shown) continue;

        var trk = tracks[i];
        tURLs[i] = trk.pURL; 
        if (trk.affinity_token) affTokens[i] = trk.affinity_token;
        trk.rg = i;
        trx[i] = '[' + SeqView.TM.trackToString(trk, null) + extraOpt + ']';
        if (trk.key == 'graph_overlay' && trk.subTracks) {
            trk.subTracks.forEach(function(uId, idx){
                for (var j = 0; j < tracks.length; j++) {
                    if (tracks[j].id != trk.legend[idx].id) continue;
                        trx[i] += '[' + SeqView.TM.trackToString(tracks[j], null)+ ', shown:false' + extraOpt + ']';
                    break;
                }
            });
        }
    }

    if (SeqView.TM.renderStat) {
        tracks.forEach(function(t) {
            console.log(t.display_name + ', render_group: ' + t.render_group + ', rg# ' + t.rg + ', order: ' + t.order);
        });
    }
    return {tracks: trx, tURLs: tURLs, affTokens: affTokens};
}


SeqView.TM.generateTracksString = function(categories){
    categories = categories || SeqView.TM.categories;
    if (!categories) return null;

    var act_tracks = SeqView.TM.collectActiveTracks(categories);
    return SeqView.TM.tracksArrayToString(act_tracks, true);
}


SeqView.TM.generateTracksStringForConfig = function(categories, filter_shown) {
    categories = categories || SeqView.TM.categories;
    if (!categories) return null;

    var str = '';
    Ext.each(categories, function(cat) {
        Ext.each(cat.subcategories, function(subcat) {
            Ext.each(subcat.tracks, function(track) {
                if (track.shown || !filter_shown) {
                    str += '[';
                    str += SeqView.TM.trackToString(track);
                    str += ',order:' + track.order;
                    str += ',shown:' + (track.shown ? 'true' : 'false');
                    str += ']';
                }
            });
        });
    });
    return str;
}


SeqView.TM.applyGroupSettings = function(categories, group_track) {
    Ext.each(categories, function(cat) {
        Ext.each(cat.subcategories, function(subcat) {
            Ext.each(subcat.tracks, function(track) {
                if (track.key === group_track.key &&
                    track.setting_group === group_track.setting_group) {
                    track.check_boxes = group_track.check_boxes;
                    track.choice_list = group_track.choice_list;
                    if (group_track.text_boxes) track.text_boxes = group_track.text_boxes;
                    if (group_track.range_controls) {
                        track.range_controls = group_track.range_controls;
                        Ext.each(track.range_controls, function(rc) {
                            var src = Ext.each(group_track.hidden_settings, function(hs) { return hs.name !== rc.name; });
                            if (src === true) return;
                            src = group_track.hidden_settings[src];
                            var idx = Ext.each(track.hidden_settings, function(hs) { return hs.name !== rc.name; });
                            if (idx !== true) track.hidden_settings[idx].value = src.value;
                            else track.hidden_settings.push(src);
                        });
                    }
                    
                }
            });
        });
    });
}

////////////////////////////////////////////////////////////////////////////////////////////////
//
SeqView.TM.getGroupSettingsString = function(group_track) {
    var str = '';
    if (group_track) {
        str += '[';
        str += SeqView.TM.trackToString(group_track, true);
        str += ']';
    }
    return str;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//
SeqView.TM.Common = {
    updateSeqViewApp: function(categories, app, group_track) {
        var group_settings;
        app = app || this.sviewApp;
        var user_data = app.m_Key;
        if (this.getEl) this.getEl().mask('Configuring...');
        var seqConfig = this.seqConfig || app.m_Config;
        categories = categories || SeqView.TM.processTracksInfo(seqConfig.TrackConfig);
        group_track = group_track || this.group_track || null;
        if (group_track) {
            SeqView.TM.applyGroupSettings(categories, group_track);
            group_settings = SeqView.TM.getGroupSettingsString(group_track);
        }
        if (!user_data || user_data.length == 0) user_data = null;
        if (this.defaultTracks) app.forEachView(function(v){ delete v.m_showAllGenes; });
        if (this.cleanCfg) this.cleanCfg();
        app.updateConfig(seqConfig, user_data);
    },

    createMainPanel: function(categories, tracks_config) {
        var removeButton = function() {
            var agp = this.tracksPanel.getActiveGridPanel();
            var win = new Ext.Window({title: 'Confirmation', modal: true, labelWidth: 1,
                constrain: true,
                app: this.sviewApp,
                layout:'form', width:460, buttonAlign: 'center',
                buttons:[
                    {text: 'Continue', handler: function() {
                        this.sviewApp.pingClick('9-2');
                        win.close();
                        this.removeUUDTracks(); }, scope: this},
                    {text: 'Cancel', handler: function() { win.close(); }}]
            });
            win.add({xtype: 'displayfield', value: 'The following track(s) will be removed from the menu after pressing "Configure":',
                style: { 'font-size': '13px' }});

            Ext.each(agp.getSelection(), function(t, idx){
                var trk = t.getData();
                if (trk.id.search(/U|R/) != 0) 
                    agp.getSelectionModel().deselect(t);
                else
                    win.add({xtype: 'displayfield', value: (idx + 1) + '. ' + trk.display_name,
                         style: { 'font-weight': 'bold', 'font-size': '13px',
                        'background-color': (idx%2 ? '' : 'gainsboro')}});
            });
            win.show();
        }

        this.categories = categories;
/*        this.activeTrackPanel = new SeqView.TM.TracksPanel({
            title: 'Active Tracks',
            iconCls: 'xsv-seq-logo',
            activeTracks: true,
            tracks_config: tracks_config,
            tracks: SeqView.TM.collectActiveTracks(categories),
            listeners: {
                deactivate: function() {
                    this.gridPanel.syncModifiedTracks(this.tracks_config);
                }
            }
        });*/
        this.tracksPanel = new SeqView.TM.CategoryTrackPanel({
            sviewApp: this.sviewApp,
            categories: this.categories
        });

        this.uploadPanel = new SeqView.TM.UploadPanel({
            cfgPanel: this,
            dirty: false
        });

        this.searchPanel = new SeqView.TM.SearchPanel({
            cfgPanel: this,
            categories: this.categories,
            dirty: false
        });

        var activeTab = this.activeTab || 0;
        var panel = new Ext.TabPanel({
            activeTab: activeTab,
            buttonAlign: 'left',
            items: [/*this.activeTrackPanel,*/ this.tracksPanel, this.searchPanel, this.uploadPanel],
            listeners: {
                beforetabchange: function(p, newtab, curtab) {
                    var isTP = newtab === this.tracksPanel;
                    Ext.getCmp('sv-delete_button_id')[isTP ? 'show' : 'hide']();
                    if (isTP && newtab.dirty) {
                        this._loadTracks(false);
                    }
                },
                scope: this
            },
            buttons: [
                {text: 'Remove track(s)', scope: this, icon: SeqView.base_url + 'images/delete.png',
                tooltip: 'Remove track(s) from NCBI data storage',
                disabled: true, id: 'sv-delete_button_id',
                handler: removeButton},
                '->',
/*                {text: 'Configure Group', scope: this,
                hidden: true, id: 'sv-config_group_button_id',
                handler: function(b) {
                    if (this.cfgWindow) {
                        if (this.sviewApp) this.group_track = b.track;
                        this.reConfigure(false);
                    }
                }},*/
                {text: 'Configure', id: 'sv-configure_button_id', scope: this, handler: function(b) {
                    this.sviewApp.cfgEvent = true;
                    this.sviewApp.pingClick('9-0');

                    if (!this.uploadPanel.uploadButt.disabled) {
                        var msgBox = Ext.MessageBox.show({title:'External Data Upload',
                            msg: 'Upload external data?',
                            buttons:Ext.MessageBox.YESNO,
                            icon:Ext.MessageBox.QUESTION,
                            fn: function(buttonId) {
                                this.reConfigure(buttonId == 'yes');
                                delete msgBox;
                            },
                            scope: this
                        });
                    }
                    else this.reConfigure(false);
                }},
                {text: 'Reset tracks', id: 'sv-load_default_button_id', scope: this, tooltip: 'Reset track list and settings',
                 handler: function(b) {
                    this.sviewApp.pingClick('9-1');
                    var msgBox = Ext.MessageBox.show({title:'Warning!',
                        msg: 'All track settings will be reset to default values. Continue?',
                        buttons: Ext.MessageBox.YESNO,
                        icon: Ext.MessageBox.WARNING,
                        fn: function(buttonId) {
                            if (buttonId == 'yes') {
                                this._loadTracks(this.defaultTracks = true);
                            }
                            delete msgBox;
                        },
                        scope: this
                    });
                }},
                {text: 'Cancel', id: 'sv-cancel_configure_button_id', scope: this, handler: function (b) {
                    this.sviewApp.pingClick('9-3');
                    this.cleanCfg();
                }}
            ]
        });
        return panel;
    },

    removeUUDTracks: function() {
        var agp = this.tracksPanel.getActiveGridPanel(),
            aStore = this.tracksPanel.activeTracksTab.gridPanel.store,
            sStore = this.searchPanel.availTracksStore;
        var sm = agp.getSelectionModel();
        var selection = sm.getSelection();
        this.removeList = selection.concat(this.removeList || []);
        Ext.each(selection, function() {
            var cfgIdx = this.getData().cfgIdx;
            aStore.removeAt(aStore.find('cfgIdx', cfgIdx));
            sStore.removeAt(sStore.find('cfgIdx', cfgIdx));
            if (sm.lastSelected == this && !sm.selectNext()) sm.selectPrevious();
            agp.store.remove(this);
        });
        if (!agp.store.getData().length) Ext.getCmp('sv-delete_button_id').disable();
    },
    // private
    _loadTracks: function(default_cfg) {
        this.getEl().mask('Loading...');
        if (default_cfg)
            this.sviewApp.m_Key = this.sviewApp.defaultConfig.m_Key;
        this.seqConfig.load({
            defaultcfg: default_cfg,
            forcereload: true,
            user_data: this.sviewApp.m_Key,
            callback: {
                success: function(config) {
                    this.getEl().unmask();
                    if (default_cfg) {
                        this.reConfigure();
                        return;
                    }
                    this.categories = SeqView.TM.processTracksInfo(config.TrackConfig);
                    this.tracksPanel.reload(this.categories);
                },
                failure: function(res) {
                    var fromcgi = Ext.decode((typeof res === 'object') ? res.responseText : res);
                    var msg = fromcgi.success === false ? fromcgi.msg : res.responseText;
                    this.getEl().unmask();
                    Ext.MessageBox.show({title: 'Configuration loading error',msg: msg,
                                         buttons: Ext.MessageBox.OK,icon:Ext.MessageBox.ERROR});
                },
                scope: this
            }
        });
    }
}

SeqView.TM.showConfigDialog = function(sview_app, activeTab) {
    sview_app.resizeIFrame(600);
    var seqconfig = new SeqView.Config(sview_app);
    seqconfig.TrackConfig = Ext.clone(sview_app.m_Config.TrackConfig);
    seqconfig.Options = Ext.clone(sview_app.m_Config.Options);
    seqconfig.SeqInfo = sview_app.m_Config.SeqInfo;
    var tm = new Ext.Window({
        title: 'Configure Page',
        constrain: true,
        app: sview_app,
//        modal: true,
        layout:'fit',
        width:800,
        height:600,
        plain: true,
        cls: 'SeqViewerApp'
    });
    var cfgPanel = Ext.create('SeqView.TM.ConfigPanel', {
       cfgWindow: tm,
       accession: sview_app.m_Config.SeqInfo.id,
       seqConfig: seqconfig,
       sviewApp: sview_app,
       activeTab: activeTab || 0
    });
    tm.add(cfgPanel);
    tm.on('close', function(p) {
        Ext.getBody().unmask();
        sview_app.resizeIFrame();
//        sview_app.m_DialogShown = false;
        sview_app.fireEvent('configuration_panel', sview_app, sview_app.cfgEvent ?  'close' : 'cancel');
        delete sview_app.cfgEvent;
    });
    sview_app.fireEvent('configuration_panel', sview_app, 'open');
//    sview_app.m_DialogShown = true;

    if (sview_app.m_iFrame) Ext.defer(tm.show, 500, tm); else tm.show();
    Ext.getBody().mask();
};


SeqView.TM.modifyTrackDetails = function(gview, track) {
    var app = gview.m_App;
    app.resizeIFrame(600);
    if (typeof track == 'number') track = app.m_Config.TrackConfig[track]; 

    var tsWindow = new Ext.Window({
        title: track.display_name,
        constrain: true,
        app: app,
        modal: true,
        width:450,
        cls: 'SeqViewerApp'
    });
    var detailPanel = new SeqView.TM.TrackDetail({track: Ext.clone(track)});

    var applySettings = function(button) {
        if (track.highlights && track.choice_list) {
            var clIdx = track.choice_list.findIndex(function(cl) { return cl.display_name == 'Rendering options'; });
            if (clIdx != -1 && track.choice_list[clIdx].curr_value != detailPanel.track.choice_list[clIdx].curr_value) delete track.highlights;
        }    
        Ext.copy(track, detailPanel.track, 'choice_list, check_boxes, range_controls, text_boxes, shown');

        if (track.category.name == 'Genes') app.m_Config.cleanGeneMode(track);

        SeqView.TM.Common.updateSeqViewApp(null, app, button.group_track);
        tsWindow.close();
    }
    tsWindow.add(detailPanel);
    var gname = track.setting_group;
    if (gname) {
        var display_limit = 17;
        if (gname.length > display_limit) 
            gname = gname.substr(0, display_limit) + "...";
    }
    var panel = tsWindow.add({buttons: [
        {text: 'Apply to Group "' + gname + '"', tooltip: 'Apply to Group "' + track.setting_group + '"',
         scope: this, handler: applySettings,
         hidden: !(gname), group_track: track},
        {text: 'Apply',  scope: this, handler: applySettings},
        {text: 'Cancel', scope: this, handler: function() {tsWindow.close();}}]});

    tsWindow.on('render', function(){
        panel.getDockedItems()[0].insert(0, [{ xtype: 'component',
            html: '<a href="' + SeqView.getHelpURL('legends/#' + track.legend_text) + '" target="_blank">Track legend</a>'}, '->']);
    });
    tsWindow.on('close', function() {
        app.resizeIFrame();
    });
    if (app.m_iFrame) Ext.defer(tsWindow.show, 500, tsWindow); else tsWindow.show();
};


Ext.define('SeqView.TM.ConfigPanel', {
    extend: 'Ext.Panel',
    initComponent : function(){
        var categories = SeqView.TM.processTracksInfo(this.seqConfig.TrackConfig);
        Ext.apply(this, {
            border: false,
            height:  600,
            layout: 'fit',
            items: this.createMainPanel(categories, this.seqConfig.TrackConfig)
        });
        this.callParent(arguments);
    },

    cleanCfg: function() {
        if (this.getEl().isMasked()) this.getEl().unmask();
        delete this.defaultTracks;
        this.cfgWindow.close();
    },

    reConfigure: function(checkDirty) {
        if (this.sviewApp) {
            if (checkDirty && !this.uploadPanel.uploadButt.disabled) {
                this.uploadPanel.submitData(this.down('#loaddatapanel').getForm(),
                                            this.cfgProcess, this);
            } else this.cfgProcess();
        } else this.cleanCfg();
    },

    cfgProcess: function() {
        var actTrx = this.sviewApp.getActiveTracks(),
            seqConfig = this.sviewApp.m_Config,
            trkCfg = this.seqConfig.TrackConfig,
            tmsTrx = this.searchPanel.checkedTracks.filter(function(t) {
                if (t.data && !t.data.hasOwnProperty('cfgIdx')) return t.data;
                if (t.hasOwnProperty('shown') && !t.shown)
                    seqConfig.TrackConfig[t.cfgIdx].shown = trkCfg[t.cfgIdx].shown = false;
            });

        this.sviewApp.addUUDTracks(tmsTrx.map(function(t) {return t.data;}));
        if (this.tracksPanel.dirty || tmsTrx.length) {
            this.getEl().mask('Loading...');
            this.seqConfig.load({
                updateUUD: true,
                forcereload: true,                
                user_data: this.sviewApp.m_Key,
                callback: {
                    success: function(config) {
                        var categories = SeqView.TM.processTracksInfo(config.TrackConfig);
                        this.updateSeqViewApp(categories);
                    },
                    failure: function(res) {
                        var fromcgi = Ext.decode((typeof res === 'object') ? res.responseText : res);
                        var msg = fromcgi.success === false ? fromcgi.msg : res.responseText;
                        delete this.sviewApp.cfgEvent;
                        Ext.MessageBox.show({title: 'Configuration loading error',msg: msg,
                                           buttons: Ext.MessageBox.OK,icon:Ext.MessageBox.ERROR});
                    },
                    scope: this
                }
            });
        } else {
            this.searchPanel.availTracksStore.getModifiedRecords().forEach(function(mt) {
                seqConfig.TrackConfig[mt.data.cfgIdx].shown = trkCfg[mt.data.cfgIdx].shown = mt.data.shown;
            });
            Ext.each(this.removeList, function(rt) { 
                var idx = seqConfig.removeUserTrack(rt.id);
                if (idx >= 0) trkCfg.splice(idx, 1);
            });
            //update cfgIdx
            if (this.removeList) {
                trkCfg.forEach(function(t, idx) { t.cfgIdx = idx; });
                seqConfig.TrackConfig.forEach(function(t, idx) { t.cfgIdx = idx; });
                this.sviewApp.forEachView(function(v) { if (v.isGraphic()) v.syncTracks(trkCfg); });
            }
/*            var agp = this.tracksPanel.getActiveGridPanel();
            if (agp.syncModifiedTracks() && this.tracksPanel.searchPanel.gridPanel == agp) {
                // if the active is Search Tracks tab and there were changes it overwrites changes made in other tabs
                var tcfg = this.tracksPanel.sviewApp.m_Config.TrackConfig;
                this.categories = SeqView.TM.processTracksInfo(tcfg);
                this.seqConfig.TrackConfig = Ext.clone(tcfg);
            }*/
            this.updateSeqViewApp(this.categories);
        }
    }
});

Ext.apply(SeqView.TM.ConfigPanel.prototype, SeqView.TM.Common);

SeqView.TM.showSeqRegion = function(title, seq) {
    var mbw = new Ext.Window({
        title: title,
        constrain: true,
        modal: true, layout:'fit',
        overflowY: 'auto',
        width: 400, height: 300,
        items:[{xtype: 'displayfield', value: seq, textalign: 'left'}],
        buttons:[{text: 'Close', handler: function() {mbw.close();}}]});
    mbw.show();
}

SeqView.TM.modifyLegendSettings = function(gview, track, idx) {
    var app = gview.m_App;
    var legends = track.legends;
    if (!legends[0].name)
        Ext.each(legends, function() {
            var split = this.label.split(': ');
            if (split.length > 2) split[1] = this.label.substr(split[0].length + 2);
            this.gStyle = split[0];
            this.name = split[1];
            this.hidden = !this.bounds;
        });
    app.resizeIFrame(300);
    var lgWindow = new Ext.Window({
        title: track.display_name,
        constrain: true,
        app: app,
        modal: true,
        width:550,
        listeners: {
            close: function() { app.resizeIFrame(); },
            show: function(){ updateSettings(); }
        }
    });
    var applySettings = function(button) {
        var doReload = false;
        var trackConfig = app.m_Config.TrackConfig;
        var subTracks = [];
        Ext.each(track.legends, function(lg) {
            if (!lg.bounds !== lg.hidden) {
                doReload = true;
                delete lg.bounds;
            }
            subTracks.push(trackConfig[lg.idx].uId + (lg.hidden ? '_hidden' : ''));
            Ext.each(trackConfig[lg.idx].hidden_settings, function() {
                if (lg[this.name] != this.value) {
                    doReload = true;
                    this.value = lg[this.name];
                }
            });
            if (lg.label.split(': ')[0] != lg.gStyle) {
                doReload = true;
                Ext.each(trackConfig[lg.idx].choice_list, function() {
                    if (this.name != 'style') return true;
                    this.curr_value = lg.gStyle;
                    return false;
                });
            }
        });
        if (doReload) {
            trackConfig[track.cfgIdx].subTracks = subTracks;
            SeqView.TM.Common.updateSeqViewApp(null, app);
        }
        lgWindow.close();
    }

    var trackStore = [];
    Ext.each(legends, function() { trackStore.push(this.name); });
    var trackPicker = new Ext.form.ComboBox({ triggerAction: 'all', width: 530, name: 'trackPicker', mode:'local',
                       store: trackStore, allowBlank: false, editable: false, value: legends[idx].name,
                       listeners: {select: updateSettings}});

    lgWindow.add(trackPicker);
    trackPicker.getTrigger('picker').setTooltip('All subtracks');

    var settingsPanel =  new Ext.Panel({
        bodyStyle: {background: 'inherit'},
        padding: '0 3 0 3'
    });
    lgWindow.add(settingsPanel);

    var colorPicker = {xtype: 'container', layout: 'hbox',  hei_ght: 22, padding: '3', // bodyStyle: {background: 'inherit'},
        items: [
           {xtype:'displayfield',  labelAlign: 'right', fieldLabel: 'Color' },
           {xtype:'component', itemId: 'colorBox', width:22, height: 22},
           {xtype:'button',  height: 22,
                menu: new Ext.menu.ColorPicker({listeners: {select: 
                    function(picker, choice) {
                        legends[idx].color = '#' + choice;
                        settingsPanel.down('#colorBox').el.applyStyles('background-color: ' + legends[idx].color + ';');
                    }, scope: this }}
                )
           }
       ]
    };
    var slidePicker = {xtype: 'container', layout: 'hbox',  height: 22, //bodyStyle: {background: 'inherit'},
        items: [
            {xtype: 'slider', itemId: 'opacity', labelAlign: 'right', fieldLabel: 'Opacity', width: '50%', increment: 5,
                listeners: {changecomplete: function(s, n, t) { this.tipText(t); }},
                tipText: function(thumb) {
                    legends[idx].opacity = thumb.value;
                    settingsPanel.down('#colorBox').el.applyStyles('opacity:' + (thumb.value/100) + ';');
                    return String(thumb.value) + '%';
            }}
       ]
    };
    var stylePicker = new Ext.form.RadioGroup({fieldLabel: 'Display style', columns: 1, labelAlign: 'right',
            items: [
                {boxLabel: 'Histogram', name: 'gStyle', inputValue: 'histogram'},
                {boxLabel: 'Line graph', name: 'gStyle', inputValue: 'line graph'}
            ],
            listeners: {change: function(sp, val) { legends[idx].gStyle = val.gStyle; }}
     });

    var trackHider = new Ext.form.Checkbox({ boxLabel: 'Hide subtrack', padding: '5',
        listeners: {change: function(o, val) { settingsPanel[(legends[idx].hidden = val) ? 'disable' : 'enable'](); }}
    });
    lgWindow.add(trackHider);

    function updateSettings(tPicker) {
        if (tPicker) {
            var rv = tPicker.getRawValue();
            for (var i = 0; i < legends.length; i++) {
                if (legends[i].name != rv) continue;
                idx = i;
                break;
            }
        }
        var legend = legends[idx];
        settingsPanel.down('#colorBox').el.applyStyles('background-color:' + legend.color + ';opacity:' + (legend.opacity/100) + ';');
        settingsPanel.down('#opacity').setValue(parseInt(legend.opacity));
        stylePicker.setValue({gStyle: legend.gStyle});
        trackHider.setValue(legend.hidden);
        // disable the checkbox checking if only one subtrack is available (not hidden)
        trackHider.setDisabled(!legend.hidden && legends.filter(function(lg){return !lg.hidden;}).length < 2);
    }

    settingsPanel.add(colorPicker, slidePicker, stylePicker);

    lgWindow.add({buttons: [
        {text: 'Apply',  scope: this, handler: applySettings},
        {text: 'Cancel', scope: this, handler: function() {lgWindow.close();}}]});

    if (app.m_iFrame) Ext.defer(lgWindow.show, 500, lgWindow); else lgWindow.show();
};

SeqView.TM.showTrackHubDialog = function(svApp) {
    var hubCnt,
        actTracks = svApp.getActiveTracks();
    var applyData = function() {
    try {
        hubCnt.getUpdatedTracks().then(function(tracks) {
            var addedTrx = [];
            for (var i = 0; i < tracks.length; i++) { // applying changes
                if (svApp.m_Config.updateOptsFromUUD(tracks[i]) < 0) addedTrx.push(tracks[i]);
            }
            svApp.addUUDTracks(addedTrx);
            svApp.m_Config.load({
                 forcereload: true,
                 updateUUD: true,
                 user_data: svApp.m_Key,
                 callback: {
                     success: function(config) {
                         SeqView.TM.Common.updateSeqViewApp(SeqView.TM.processTracksInfo(config.TrackConfig), svApp);
                     },
                     failure: function(res) {
                         var fromcgi = Ext.decode((typeof res === 'object') ? res.responseText : res);
                         var msg = fromcgi.success === false ? fromcgi.msg : res.responseText;
                         delete svApp.cfgEvent;
                         Ext.MessageBox.show({title: 'Configuration loading error', msg: msg,
                             buttons: Ext.MessageBox.OK,icon:Ext.MessageBox.ERROR});
                     }
                }
            });
        });
    } catch(e) { console.log(e); tm.close();}
    }
    GB.HubContent.setApp(svApp);
    if (GB.Panel != null) {
        hubCnt = GB.HubContent.showHubs({assm: svApp.m_AssmContext, commit: applyData, tracks: actTracks}).getContent();
        return;
    }
    hubCnt = new GB.HubContent();
    var id = 'hubCnt' + svApp.m_Idx;
    var tm = new Ext.Window({
        title: 'Configure Track Hubs',
        constrain: true,
        modal: true,
        width:800,
        height:600,
        buttons: [
            {text: 'Configure', handler: function() { applyData(); tm.close();}},
            {text: 'Cancel', handler: function() { hubCnt.revert(); tm.close()}}
        ],
        cls: 'SeqViewerApp'
    });
    svApp.resizeIFrame(tm.height);
    tm.panel = tm.add({
        html: {tag: 'div', id: id},
        layout: 'fit',
        height: '100%',
        region: 'center'
    });
    tm.on('resize', function() {
        this.panel.updateLayout();
        this.panel.setHeight(this.body.getHeight());
        hubCnt.resizeContentView({height: this.body.getHeight() - 4});
    });
    tm.on('close', function(p) {
        svApp.resizeIFrame();
        svApp.fireEvent('trackhub_panel', svApp, svApp.cfgEvent ?  'close' : 'cancel');
        delete svApp.cfgEvent;
    });
    tm.on('show', function(p) {
        document.getElementById(id).appendChild(hubCnt.view);
        this.updateLayout();
        hubCnt.loadContent({assm: svApp.m_AssmContext, tracks: actTracks});
    });

    svApp.fireEvent('trackhub_panel', svApp, 'open');
    if (svApp.m_iFrame) Ext.defer(tm.show, 500, tm); else tm.show();
};


/*  $Id: tm_grids.js 46495 2021-06-07 19:56:47Z borodine $
*/

Ext.define('SeqView.TrackList', {
    extend: 'Ext.data.Model',
    proxy: {type: 'memory'}, //enablePaging: true},
    fields: [
        'display_name',
        'help',
        { name: 'order', type: 'int' },
        { name: 'shown', type: 'bool' },
        { name: 'subcat', mapping: function(data)
            {return data.subcat || data.category.display_name;} }
    ]
});

Ext.define('SeqView.TM.GridBase', {extend: 'Ext.grid.GridPanel',
    border: false,
    itemId: 'gridbase',
    multiSelect: true,
    overflowY: 'auto',
    region: 'center',
    initComponent: function() {
        if (this.subcategories && this.subcategories.length > 1) {
        // groupping tracks panel
            var tracks = [];
            Ext.each(this.subcategories, function(subcat) {
                Ext.each(subcat.tracks, function(track) {
                    track.subcat = (1000000 + subcat.order).toString() + subcat.display_name;
                    if (SeqView.debugMode) track.subcat += ' ' + subcat.order;
                    tracks.push(track);
                });
            });
            this.tracks = tracks;
            this.groupGrid = true;
        }
        if (this.searchTab) this.tracks.forEach(function(t) { delete t.subcat; }); 
        this.store = this.store || Ext.create('Ext.data.Store', {
            autoDestroy: true,
            model: 'SeqView.TrackList',
            data: this.tracks
        });
        if (this.groupGrid) {
            this.store.setGroupField ('subcat');
            this.features = [{
                ftype: 'grouping',
                groupHeaderTpl: ['Category: ', '{name:this.subcatName} ({[values.rows.length]} Item{[values.rows.length > 1 ? "s" : ""]})',
                    { subcatName: function(name) { return name.substr(7); } }
                ],
                hideGroupedHeader: true,
                enableGroupingMenu: false
            }];
        }
        
        var sortable = this.hasOwnProperty('sortable') ? this.sortable : false;

        if (!this.columns) {
            this.columns = [
                {header: 'Active', width: 45, xtype: 'checkcolumn', sortable: sortable, dataIndex: 'shown', menuDisabled: true,
                    listeners: {scope: this, checkchange: this.checkChange},
                    editor: { xtype: 'checkbox', cls: 'x-grid-checkheader-editor'}},
                {header: "Track name", sortable: sortable, dataIndex: 'display_name', flex: true, editable: false, menuDisabled: true}];
                if (this.searchTab) {
                    this.columns.unshift({header: "Category", width: 70, sortable: true, dataIndex: 'subcat', editable: false, menuDisabled: true});
                    this.columns.push({header: "Description", dataIndex: 'help', flex: true, editable: false, menuDisabled: true});
                }
        }
        if (this.groupGrid) this.columns.push({header: 'Category', dataIndex: 'subcat', sortable: sortable, width: 70});
        this.callParent(arguments);
    },
    
    checkChange: function(col, num, val, trk) {
        var sel = this.getSelection();
        if (sel.find(function(t){ return t == trk; }))
            sel.forEach(function(t){ t.set('shown', val); });
    },
    
    syncModifiedTracks: function(tracks_config) {
        var mr = this.getStore().getModifiedRecords();
        tracks_config = tracks_config || this.tracks;
        for (var i = 0; i < mr.length; i++) {
            var data = mr[i].data;
            Ext.each(tracks_config, function(tr) {
                if (tr.id == data.id) { tr.shown = data.shown; return false; }
            });
        }
        return mr.length;
    }
});


SeqView.TM.createSearchTracksTab = function(title, tracks_config, str) {
    var trackFilter = function(rec) {
        if (!this.baseParams.regexp) return true;
        try {
            if (rec.data.help.search(this.baseParams.regexp) >= 0) return true;
            if (rec.data.display_name.search(this.baseParams.regexp) >= 0) return true;
        } catch(err) { console.log(err.message); }        
        return false;
    };
    var pageSize = 10; 
    var panel = new SeqView.TM.GridBase({
        title: (title) ? title : null,
        border: false,
        searchTab: true,
        sortable: true,
        stripeRows: true,
        tracks: tracks_config,
        viewConfig: {forceFit:true, deferEmptyText:false, emptyText:'<div align="center">No Search Results To Display</div>'}/*,
        bbar: new Ext.PagingToolbar({
            store: store,
            displayInfo: true, displayMsg: 'Displaying Search Results {0} - {1} of {2}', emptyMsg: "No Search Results To Display"
        })*/
    });
    var store = panel.store;
    store.baseParams = {query: str || '', limit: pageSize, filterBy: trackFilter, regexp: str ? globStringToRegex(str) : ''};
    if (str) store.filterBy(trackFilter, store);
    if (title) panel.iconCls = (store.getCount()) ? 'xsv-search-results' : '';
    return panel;
};


/*  $Id: sviewapp.js 47512 2023-06-13 20:59:13Z asztalos $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Authors:  Vlad Lebedev, Maxim Didenko, Victor Joukov
 *
 * File Description:
 *
 */

SeqView.fireEvent = function(eName, elemID){
    var elem = document.getElementById(elemID);
    if (document.createEvent) {
        var e = document.createEvent('HTMLEvents');
        e.initEvent(eName, false, false);
        elem.dispatchEvent(e);
    } else elem.fireEvent(eName);
};


SeqView.Config = (function() {

    function constructor(app) {
        this.m_App = app;
    }
    return constructor;

}) ();

SeqView.Config.prototype = {
    checkWarnings: function() {
        var wrn = '';
        Ext.each(this.SeqInfo.warning_messages || [], function() { wrn += this + '<br>'; });
        this.SeqInfo.warnTT = wrn;
        this.SeqInfo.icon = SeqView.ExtIconLoading.slice(5, -18) + (wrn ? 'shared/warning.gif' : 'window/toast/icon16_info.png');
        return wrn;
    },

    trackInfo2CfgArray: function(tInfo, cfgIdx) {
        cfgIdx = cfgIdx || 0;
        var trx = [];
        tInfo = tInfo || this.trackInfo;
        tInfo.forEach(function(t) {
            t = t.track || t;
            var tt = {id: t.id, key: t.key, shown: t.shown != 0 || t.shown != false};
            Ext.each(['annots', 'batch', 'big_data_index', 'comments', 'category', 'check_boxes', 'choice_list', 'data_key',
                    'data_provider',
                    'dbname', 'display_name', 'extra_info', 'filter', 'help', 'hidden_settings', 'highlights', 'highlights_color',
                    'hub_id', 'hub_name', 'hub_url', 'is_private', 'legend', 'legend_text', 'name', 'order', 'range_controls',
//                    'remote_data',
                    'rmt_mapped_id', 'seq_data_key', 'setting_group', 'stored_scale', 'subcategory', 'subkey', 'subTracks', 'text_boxes'
            ], function(p) { var v = t[p]; if (v) tt[p] = v; });
            Ext.each(tt.check_boxes, function(cb) { cb.value = (typeof cb.value == 'boolean' ? cb.value : cb.value == 'true'); });
            tt.cfgIdx = cfgIdx++;
            if (t.shownBool) tt.tmsSearch = true;
            trx.push(tt);
        });
        return trx;
    },

    set_pURL_uId: function(firstIdx, last) {
        var trx = this.TrackConfig,
            last = last || trx.length,
            baseURL = this.m_App.m_CGIs.baseURL;
        for (var i = firstIdx; i < last; i++) {
            var t = trx[i];
            t.pURL = baseURL + this.setProviderAndAffinity(t) + '.cgi';
            Ext.each(t.subTracks, function(tId) {
                (trx.find(function(tc) { return tId == tc.id; }) || {}).uId = tId;
            });
        }
    },

    updateOptsFromUUD: function(uudTrack) {
        var cfgIdx = this.TrackConfig.findIndex(function(t) { return t.id == uudTrack.id; });
        if (cfgIdx >= 0) {
            if (uudTrack.active === null) {
                this.TrackConfig.splice(cfgIdx, 1);
                return 0;
            }
            var tt = this.TrackConfig[cfgIdx];
            tt.shown = uudTrack.shown != 0 || uudTrack.shown != false;
            Ext.each(['check_boxes', 'choice_list', 'hidden_settings', 'highlights', 'highlights_color', 'range_controls', 'stored_scale', 'text_boxes'
            ], function(p) { var v = uudTrack[p]; if (v) tt[p] = v; });
            Ext.each(tt.check_boxes, function(cb) { cb.value = (typeof cb.value == 'boolean' ? cb.value : cb.value == 'true'); });
        }
        return uudTrack.active === null ? 0 : cfgIdx;
    },

    load: function(cfg) {
        var cleanPanel = function(panel) {
            if (panel.tmpPanel) {
                delete panel.tmpPanel;
                panel.removeAll();
            }
        }
        this.m_App.m_timeStamps.prior2seqconfig = new Date().getTime();
        if (window.timeStamp) {
            this.m_App.ping({'SV_start_to_seqconfig_time': this.m_App.m_timeStamps.prior2seqconfig - timeStamp,
                          'sv-event':'initialization'});
            delete timeStamp;
        }

        if (this.m_App.m_Panel.tmpPanel) this.m_App.m_Panel.tmpPanel.mask('Loading...');
        try {
            // getting path to ExtJS loading.gif
            SeqView.ExtIconLoading = getComputedStyle(document.querySelector('.x-mask-msg-text'))['background-image'];
        } catch(e) {};
        if (!cfg.callback) {
            console.error('Undefined callback object in ', cfg);
            return;
        }
        Ext.applyIf(this.tmpCallback = cfg.callback, {scope: this, success: Ext.emptyFn, failure: Ext.emptyFn});

        var url = cfg.config_url || this.m_App.m_CGIs.Config;
        var seqid = cfg.seqid || this.m_App.GI;

        var params = { id: seqid };
        if (this.m_App.defaultConfig.nocache) params.nocache = this.m_App.defaultConfig.nocache;
        // Save initial parameters for Load Defaults, or restore them from saved
        // All other relevant parameters seem to be untouched
        if (cfg.defaultcfg) {
            // Restore
            this.m_App.m_TracksFromURL = this.m_App.defaultConfig.m_TracksFromURL;
            this.m_App.m_ViewContent   = this.m_App.defaultConfig.m_ViewContent;
        } 
        Ext.each([
            ['m_TracksFromURL', 'tracks'],
            ['m_ViewTheme', 'theme'],
            ['m_ViewContent','content'],
            ['m_ItemID', 'itemID'],
            ['m_AppContext', 'app_context'],
            ['m_AssmContext', 'assm_context'],
            ['m_noGuessAssm', 'noguess_assm'],
            ['m_SRZ', 'srz'],
            ['m_BamPath', 'bam_path'],
            ['m_AppName', 'appname'],
            ['m_DepthLimit', 'depthlimit']],
            function(pair) {
                if (this.m_App[pair[0]]) params[pair[1]] = this.m_App[pair[0]];
            },
            this
        );
        var cookie_name = SeqView.Cookies.UserTracksCookieName;
        if (this.m_App.m_AppName && this.m_App.m_AppName.length > 0 && cookie_name === SeqView.Cookies.UserTracksCookieNameBase) {
            cookie_name += '-' + this.m_App.m_AppName;
        }

        if (cfg.resettracks || (cfg.defaultcfg && !this.m_App.m_TracksFromURL)) {
            delete params.tracks;
            params.clear = cookie_name;
        }
        else if (!params.tracks) {
            params.tracks = !this.TrackConfig ? '[amend]'
                            : SeqView.TM.generateTracksStringForConfig(SeqView.TM.processTracksInfo(this.TrackConfig));//, true);
        }

        if (cfg.defaultcfg !== true){
            Ext.apply(params, this.visualOptionsUrl());
            var user_config_key = cfg.user_config_key || SeqView.UserTracks.get(cookie_name, null);
            if (user_config_key) params.userconfigkey = user_config_key;
        } else {
            // we need to pass app specific data in case we have one (ISCA, A2A browsers);
            var app_data = SeqView.SessionData.get(SeqView.Cookies.AppDataCookieName, null);
            if (app_data) params.key = app_data;
            if (this.Options && this.Options.curr_content) {
                params.content = this.Options.curr_content;
            }
            params.nocache = 1;
        }

        var viewer_context = cfg.viewer_context || this.m_App.m_ViewerContext;
        if (viewer_context) params.viewer_context = viewer_context;
        var user_data = cfg.user_data || this.m_App.m_Key;
        if (user_data) params.key = user_data;
        
        var paramsFast = null;

        if (this.SeqInfo && this.trackInfo) {
            this.TrackConfig = this.trackInfo2CfgArray();
            this.set_pURL_uId(0);
            delete this.trackInfo;
            this.SeqInfo.views = ['panorama', 'graphical:0-' + (this.SeqInfo.length - 1)];
//            if ((this.m_App.m_Panel.down('#btnGeneMode') || {}).pressed) this.changeGeneMode(true);
            cleanPanel(this.m_App.m_Panel);
            this.tmpCallback.success.call(this.tmpCallback.scope, this);
            return;
        }
        if (cfg.updateUUD) {
            if (this.parseUUDTracks(params)) {
                this.tmpCallback.success.call(this.tmpCallback.scope, this);
                return;
            }
            params.fastconfig = 1;
        } else {
            if (this.m_App.uudTracks) this.parseUUDTracks(params, false);
        }

        if (!this.SeqInfo) {
            this.tmpCallback.completeLoading = this.m_App.completeCfgLoading;
            if (this.m_App.m_TracksFromURL && (this.TrackConfig = this.m_App.parseTracksFromURL())) {
                //request sequence info only
                paramsFast = {notrackconfig: 1, userconfigkey: user_config_key, id: params.id, assm_context: params.assm_context }
            } else {
                //request sequence info with active track list  
                paramsFast = Ext.apply({fastconfig: 1}, params);
            }
            cfg.forcereload = true;
        }

        // We only use content param from URL ones for the first load.
        // After that we will get it for m_Config.Options.curr_content
        delete this.m_App.m_ViewContent;
        // Same. Get track configuration from TrackConfig
        delete this.m_App.m_TracksFromURL;

        if (this.TrackConfig && !cfg.forcereload) params.notrackconfig = 1;

        var req_num = this.m_App.m_ReqNum || 0;
        this.m_App.m_ReqNum = ++req_num;

        var delay = 0;
        var rqOptions;
        var checkJobStatus = function(data, text, res) {
            if (this.m_App.m_ReqNum > req_num) return;
            var from_cgi = SeqView.decode(data);
            if (typeof from_cgi !== 'object') {
                cleanPanel(this.m_App.m_Panel);
                this.tmpCallback.failure.apply(this.tmpCallback.scope, arguments);
                return;
            }
            if (from_cgi.job_status) {
                var st = from_cgi.job_status;
                if (st == 'submitted' || st == 'running' || st == 'pending') {
                    rqOptions.url = url + '?job_key=' + from_cgi.job_id;
                    Ext.defer(SeqView.App.simpleAjaxRequest,
                        Math.min(15000, Math.max(2000, 1000 * delay++)),
                        this, [rqOptions]);
                    return;
                }
            }
            cleanPanel(this.m_App.m_Panel);
            this.m_App.watchdogStop();
            if (from_cgi.success === false || from_cgi.job_status == 'failed') {
                if (!this.m_App.m_NoDataCookie && params.key ){
                    var cookie_name = SeqView.Cookies.UserDataCookieName;
                    if (this.m_App.m_AppName && this.m_App.m_AppName.length > 0 && cookie_name === SeqView.Cookies.UserDataCookieNameBase) {
                        cookie_name += '-' + this.m_App.m_AppName;
                    }
                    SeqView.SessionData.set(cookie_name, params.key);
                }
                this.tmpCallback.failure.apply(this.tmpCallback.scope, arguments);
            } else {
                var baseURL = this.m_App.m_CGIs.baseURL;
                if (!this.SeqInfo) {
                    this.SeqInfo = from_cgi.SeqInfo;
                    TMS.Direct.getContigInfo(this.SeqInfo.assm_context, [this.SeqInfo.id], {map: 1})
                        .then(function(dict) { this.m_App.m_ViewParams.mappedGDV = dict[this.SeqInfo.id]; }.bind(this)).catch(function(e) {});
                }
                this.m_App.defaultAffinityToken = this.SeqInfo.top_level ? ('toplev_' + this.SeqInfo.id) : ''; 
                if ((paramsFast || {}).notrackconfig) {
                    paramsFast = null;
                    this.TrackConfig.forEach(function(t) {
                        t.id = t.id || '';
                        var rmt = !t.id.search(/R|H/);
                        if (t.id.indexOf('STD') == 0 && t.annouts) {
                            rmt = t.annotes[0].search(/(https?|ftp):\/\//) == 0;
                        }
                        t.pURL = baseURL + this.setProviderAndAffinity(t, rmt) + '.cgi';
                        if (t.key == 'sequence_track') this.m_App.defaultAffinityToken = t.affinity_token;
                    }, this);
                    this.tmpCallback.success.call(this.tmpCallback.scope, this);
                    requestConfigData();
                    return;
                }
/*              error: function(data, text, res) {
                    if ((this.SeqInfo || {}).requested) delete this.SeqInfo;*/

                Ext.apply(from_cgi.Options, window.JSON.parse(localStorage.getItem('NCBI/SV/Preferences')));
                if (this.fastconfig) 
                    this.TrackConfig.forEach(function(ft) {
                        var idx = from_cgi.TrackConfig.findIndex(function(t) { return t.id == ft.id; });
                        if (idx < 0) return;
                        ft.order = from_cgi.TrackConfig[idx].order;
                        from_cgi.TrackConfig[idx] = ft;
                    }, this);
                delete this.fastconfig;
                Ext.each(from_cgi.TrackConfig, function(t, idx) {
                    t.pURL = baseURL + (t.data_provider || 'seqgraphic') + '.cgi';
                    if (t.key == "graph_overlay") {
                        Ext.each(t.subTracks, function(tId) {
                            if (Ext.each(from_cgi.TrackConfig, function(tc) {
                            return tId != tc.uId; }))
                                Ext.each(from_cgi.TrackConfig, function(tc) {
                                    if (tId != tc.id) return;
                                    tc.uId = tId;
                                    return false;
                                });
                        });
                    }
                    t.cfgIdx = idx;
                });
                Ext.apply(this, from_cgi);
                if (paramsFast) {
                    this.fastconfig = paramsFast.fastconfig;
                    params.tracks = this.m_App.getActTracks();
                    paramsFast = null;
                    this.tmpCallback.success.call(this.tmpCallback.scope, this);
                    requestConfigData();
                    return;
                }
                this.m_App.adjustViewOptions();
                (this.tmpCallback.completeLoading || this.tmpCallback.success).call(this.tmpCallback.scope, this);
                delete this.tmpCallback;
            }
        }

        var requestConfigData = function() {
            delay = 0;
            rqOptions = {url: url, context: this, data: paramsFast || params, success: checkJobStatus, error: checkJobStatus};

            this.m_App.watchdogStart(rqOptions.url, '', rqOptions.data);
            this.m_App.AjaxRequest(rqOptions);

        }.bind(this);
        requestConfigData();
        
    },

    visualOptionsUrl: function() {
        var params = {};
        if (this.Options) {
            Ext.each(['label','color', 'decor', 'spacing', 'content'], function(p) {
                var pattr = 'curr_'+p;
                if (this.Options[pattr] && this.Options[pattr].length > 0)
                    params[p] = this.Options[pattr];
            },this);
        }
        return params;
    },

    parseUUDTracks: function(params, toCfg) {
        toCfg = toCfg != undefined ? toCfg : true;
        var order = (this.TrackConfig || {}).length || 1000;
        this.m_App.uudTracks.forEach(function(t) {
            toCfg &= (t.track !== undefined && t.track.category !== undefined);
        });
        if (!toCfg) {
            var id_regex = /(\[|,){1}id:\w+(\]|,){1}/;
            this.m_App.uudTracks.forEach(function(track) {
                var trkParams = '';
                track.order = track.order || order++;
                Ext.each(['subtracks', 'shown', 'order'], function(prop) {
                    if (track[prop] != undefined) trkParams += ',' + prop + ':' + track[prop];
                });
                if (track.display_settings) trkParams += ',' + track.display_settings;
                var found = false;
                var name = track.name || track.GetName();
                var key = track.track_type || track.GetAttr('track_type');
                if (name !== "" && key !== "") {
                    name = 'name:' + name;
                    key = 'key:' + key;
                    var tr = params.tracks.split(/\[(.*?)\]/).filter(Boolean);
                    for (var i = 0; i < tr.length; i++) {
                        var p = tr[i].split(',');
                        if (p.indexOf(name) != -1 && p.indexOf(key) != -1) {
                            if (!id_regex.test(tr[i])) tr[i] += ',id:' + track.id;
                            tr[i] += trkParams;
                            found = true;
                        }
                    }
                }
                if (found) params.tracks = '[' + tr.join('][') + ']';
                else params.tracks +=  '[id:' + track.id + trkParams + ']';
            });
        } else {
            var uudTrx = [],
                firstIdx = this.TrackConfig.length; 
            this.m_App.uudTracks.forEach(function(trk) {
                var trk = trk.track;
                trk.order = trk.order || order++;
                var cfgIdx = this.TrackConfig.findIndex(function(t) { return t.id == trk.id; });
                if (cfgIdx < 0) uudTrx.push(trk);
                else {
                    trk.order = this.TrackConfig[cfgIdx].order;
                    this.TrackConfig[cfgIdx] = this.trackInfo2CfgArray([trk], cfgIdx)[0];
                    this.set_pURL_uId(cfgIdx, cfgIdx + 1);
                }
            }, this);
            this.TrackConfig = this.TrackConfig.concat(this.trackInfo2CfgArray(uudTrx, firstIdx));
            if (uudTrx.length) this.set_pURL_uId(firstIdx);
        }
        delete this.m_App.uudTracks;
        return toCfg; 
    },

    cleanGeneMode: function(t) {
        var opts = t.choice_list[0],
            vs = opts.values,
            name = opts.prev_value;
        delete opts.prev_value;
/*        for (var i = vs.length; i--;) {
            if (vs[i].display_name.charAt(0) != '*') continue;
            vs[i].display_name = vs[i].display_name.replace(/^\* /,'');
            name = vs[i].name;
        }*/
        return name;
    },

    changeGeneMode: function(showAll) {
        Ext.each(this.TrackConfig, function(t) {
            if (t.category.name != 'Genes' || !t.choice_list) return;
            var cur = t.choice_list[0];
            if (showAll) {
                if (Ext.each(cur.values, function(v) { return v.name != 'ShowAll' }) === true) return;
                Ext.each(cur.values, function(v) {
                    if (v.name !== cur.curr_value) return;
/*                    if (cur.prev_value && cur.prev_value == 'MergeAll') {
                        console.log('Overwrite');
                    }*/
                    cur.prev_value = v.name;
//                    v.display_name = '* ' + v.display_name;
                    return false;
                });
                cur.curr_value = 'ShowAll';
            } else cur.curr_value = this.cleanGeneMode(t) || (cur.curr_value == 'ShowAll' ? 'MergeAll' : cur.curr_value);
        }, this);
        this.m_App.fireEvent('configuration_changed', this.m_App);
    },


    setProviderAndAffinity: function(track, rmt) {
        var isRemote = function(t) {
            return getProvider(t).search(/_rmt$/) > 0;
        }
        var getProvider = function(t) {
            return t.data_provider || setProvider(t);
        }
        var setProvider = function(t) {
            var rmt = (rmt === undefined) ? t.id.search(/R|H/) == 0 : rmt;
            return t.data_provider = (t.data_provider || 'seqgraphic').replace(/_rmt$/, '') + (rmt ? '_rmt' : '');
        }
        
        if (track.affinity_token) return getProvider(track);

        if (isRemote(track)) {
	       track.affinity_token = 'fnamecs_' + NCBIGBUtils.MD5(track.annots[0]);
	       return track.data_provider;
        }
        
        // graph overlay/composite tracks get affinity token of one of its component tracks:
        // prioritization: if there is at least one that is remote, it will get its data provider and affinity
        if (track.subTracks) {
            for (var i = 0; i < track.subTracks.length; i++) {
                var id = track.subTracks[i];
                var sbt = this.TrackConfig.find(function(t) { return t.id == id; });
                if (!sbt) continue;
                track.data_provider = this.setProviderAndAffinity(sbt);
                track.affinity_token = sbt.affinity_token;
                if (isRemote(track)) break;
            }
        }
        if (!track.affinity_token) {
            var seq = this.SeqInfo;
            if (seq.top_level || (seq.placements && seq.placements.length))
                track.affinity_token = 'toplev_' + (seq.top_level ? seq.id : seq.placements[0].accession);
            else // type is like "protein", "mRNA", "genomic", "RefSeqGene" etc returned in SeqInfo.type
                track.affinity_token = 'type_' +  seq.type;
        }
        return track.data_provider;
	},

    removeUserTrack: function(id) {
        var idx = this.TrackConfig.findIndex(function(t) { return t.id == id; });
        if (idx >= 0) {
            TMS.Direct.removeTrack(this.TrackConfig[idx]).catch(console.error);
            this.TrackConfig.splice(idx, 1);
        }
        return idx;
    }
};

/********************************************************************/
//////////////////////////////////////////////////////////////////////
// SeqView.Api
/********************************************************************/
SeqView.Api = (function() {

    function constructor(app) {
        this.m_App = app;
    }

    return constructor;

}) ();

SeqView.Api.prototype = {
    getAllDisplayOptions: function() { return this.m_App.m_Config.DisplayOptions; },
    getTracks: function(all) { return this.m_App.getTrackObjects(all); },
    setTracks: function(tracks) {}
};

/********************************************************************/
//////////////////////////////////////////////////////////////////////
// SeqView.App
/********************************************************************/

Ext.define('SeqView.App', {
    extend: 'Ext.util.Observable',

    constructor: function(div_id) {
        var sm_ResizeWatcher = null;
        if (this.m_domDiv = Ext.get(div_id))
            this.m_domDiv.addCls("SeqViewerApp");
        else
            throw "A div element containing SeqViewer should have an id attribute";

        var onWindowResize = function(w, h) {
            var rw = w;
            var rh = h;
            if (sm_ResizeWatcher) {
                clearTimeout(sm_ResizeWatcher);
                sm_ResizeWatcher = null;
            }
            sm_ResizeWatcher = Ext.defer(function(w, h) {
                SeqView.m_Apps.forEach(function(app) { app.doWindowResize(w, h); });
            }, 500,this,[rw,rh]);
        }

        Ext.on('resize', onWindowResize);

        this.m_DivId = div_id;
        this.m_DivCtrlId = this.m_DivId + 'Controls';
        this.m_DivJSId = this.m_DivId + 'JS';
        this.m_DivTitle = this.m_DivId + 'Title';
        this.m_DivTitleID = this.m_DivId + 'TitleID';
        this.m_Views = [];
        this.m_timeStamps = {},
        this.m_MarkersInfo = null;
        this.m_TextView = null;
        this.m_InitialLoading = false;
        this.m_DialogShown = false;
        this.m_Idx = SeqView.m_Apps.push(this) - 1;
        this.callParent(arguments);
    },

    pingClick: function(area, svevent) {
            var obj = {jsevent: 'click', 'sv-area': area};
            if (svevent) obj['sv-event'] = svevent;
            this.ping(obj);
    },
    
    tmp_pingArgs: [],
    ping: function(a) { // saving pings while NCBI instrumented page is loading
        var args = this.tmp_pingArgs;
        if (typeof ncbi !== 'undefined' && ncbi.sg && ncbi.sg.ping) {
            this.ping = function(obj) {
                obj.ncbi_app = 'sviewer-js';
                obj['sv-appname'] = this.m_AppName;
                ncbi.sg.ping(obj);
            };
            if (SeqView.doInitPing) {
                this.ping({'sv-browser': Ext.browser.name,
                    'sv-browser-version': Ext.browser.version.major, 'sv-os': Ext.os.name});
                delete SeqView.doInitPing;
            }
            while (args.length) this.ping(args.shift());
            delete this.tmp_pingArgs;
            this.ping(a);
        } else {
            if (args.push(a) > 5) console.log('DROPPED PING: ' + JSON.stringify(args.shift()));
        }
    },
    
    getApi: function() { return new SeqView.Api(this); },
    getAppContext: function() { return this.m_AppContext; },
    getWidth: function() { return this.m_domDiv.getWidth() - (this.m_Embedded ?  0 : 2); },
    getDataLength: function() { return this.m_SeqLength; },

    showPrintPageDlg: function(app_idx) {
        var app = this.findAppByIndex(app_idx);
        if (app) { app.showPrintPageDlg(); }
    },

    doWindowResize: function(w,h) {
        this.updateLayout();
        if (this.m_BrowWidth != w) {
            this.reloadAllViews({resize: true});
            this.m_BrowWidth = w;
        }

    },

    reload: function(appParams, keep_config) {
        if (this.m_InitialLoading){
            if (this.m_deferredReload) clearTimeout(this.m_deferredReload);
            this.m_deferredReload = Ext.defer(this.reload, 100, this, [appParams, keep_config]);
            return;
        }
        delete this.m_deferredReload;
        if (this.m_AjaxTrans) {
            // kill all active ajax requests
/*            try { //doesn't work while we use jQuery AJAX
                Ext.each(this.m_AjaxTrans,function(tid) {
                    Ext.Ajax.abort(tid);
                });
            } catch(e) {}*/
            delete this.m_AjaxTrans;
        }
        this.m_Panorama = null;
        this.forEachView(function(v) { v.destroy(); });
        this.m_Views = [];

        var div = Ext.get(this.m_DivJSId);
        if (div) div.remove();
        div = Ext.get(this.m_DivCtrlId);
        if (div) div.remove();

        if (keep_config !== true)
            delete this.m_Config;
        this.load(appParams);
    },

//////////////////////////////////////////////////////////////////////////
// initLinks: initialize GI and data independent links, all other links
// are initialized in loadAccession
    initLinks: function() {
        var baseURL = this.m_baseURL || SeqView.base_url;
        this.m_CGIs = {
            baseURL:     baseURL,
            UUD_svdl: UUD.GetUUDCgiUrl('sv-dl'),
            UUD_dlmeta: UUD.GetUUDCgiUrl('dl-meta'),
//            TMS_tracksearch: TMS.GetCgiUrl('tracksearch'),
            FeatSearch: baseURL + 'sv_search.cgi',
            ObjInfo:    baseURL + 'seqgraphic.cgi',
            NetCache:   baseURL + 'ncfetch.cgi',
            Feedback:   baseURL + 'feedback.cgi',
            Config:     baseURL + 'seqconfig.cgi',
            Panorama:   baseURL + 'seqgraphic.cgi',
            Graphic:    baseURL + 'seqgraphic.cgi',
            GraphicRmt: baseURL + 'seqgraphic_rmt.cgi',
            Alignment:  baseURL + '../msaviewer/cgi/alnmulti.cgi',
            ObjCoords:  baseURL + 'seqgraphic.cgi',
            Sequence:   baseURL + 'seqfeat.cgi',
            SequenceSave: baseURL + 'sequence.cgi',
            Watchdog:   baseURL + 'sv_watchdog.cgi',
            SearchMru:  baseURL + 'sv_mru.cgi'
        };
        if (this.m_hiddenOptions.search('use_jsonp') >= 0)SeqView.jsonType = 'JSONP';
        else {
            SeqView.jsonType = 'JSON';
            if (Ext.isIE) {
                SeqView.App.simpleAjaxRequest({url: this.m_CGIs.Feedback,
                    success: function(data, text) { SeqView.jsonType = 'JSON'; },
                    error: function(data, text) { SeqView.jsonType = 'JSONP'; }
                });
            }
        }
        
        if (!this.m_AppName || this.m_AppName == 'no_appname') {
            if (!SeqView.origHostName && !this.m_Embedded && SeqView.standalone) this.m_AppName = 'ncbi_sviewer';
            else {
                var newAppName = document.location.hostname;
                if (SeqView.origHostName || (newAppName.search(/www|qa|test|blast/) != 0)) {
                    var msg = 'Sequence Viewer was initialized with parameter <b>appname=' + this.m_AppName + '</b><br>';
                    Ext.Msg.show({title: 'Warning',
                        msg: msg + 'It will be substituted with the hostname: <b>' + newAppName + '</b>',
                        icon: Ext.Msg.WARNING,
                        buttons: Ext.Msg.OK
                    });
                }
                this.m_AppName = newAppName;
                this.m_appWarning = 'Uninitialized <b>appname</b> is set to <b>' + newAppName + '</b><br>';
           }
        }
        if (TMS.setContext) TMS.setContext(this.m_AppContext);
    },

    loadSequence: function(seqInfo, trackInfo, appParams) {
        if (this.m_MarkersInfo) this.m_MarkersInfo.forEachMarker(function(m){
            m.deleteMarker();
        });
        this.load({SeqInfo: seqInfo, trackInfo: trackInfo, appParams: appParams});
    },

    load: function(cfg, maxDelay) {
        maxDelay = maxDelay || 1000;
        var interval = 3;
        var app = this;
//        fetch('data/display_options.json', {mode:'cors'}).then(function(r) { return r.json(); }).then(function(data) {  console.log(data); });
        setTimeout(function() { // Delay to allow jQuery to get ready
            maxDelay -= interval;
            if (typeof TMS === 'undefined' && maxDelay > 0)
                app.load(cfg, maxDelay);
            else {
                var loadSeq = false,
                    appParams = cfg;
                if (cfg) {
                    loadSeq = cfg.SeqInfo != undefined;
                    appParams = cfg.appParams || cfg;
                }
                cfg = cfg || {};

                if (!loadSeq) {
                    if (app.m_Config) {
                        console.warn('SeqView.App.load() can be invoked only once per app instance. Use SeqView.App.reload()');
                        return;
                    }   
                }
                app.clearParams(loadSeq);
                if (loadSeq) app.GI = cfg.SeqInfo.id;
                app.parseParams(appParams, document.location.href.split("?")[1]);
                app.initLinks();
                app.m_Config = Ext.apply(new SeqView.Config(app), cfg);
                
                if (loadSeq) {
                    if (app.m_AppName != 'dataviewer')
                        TMS.Direct.getContigInfo(cfg.SeqInfo.assm_context, [app.GI], {map: 1})
                            .then(function(dict) { app.m_ViewParams.mappedGDV = dict[app.GI]; }.bind(this)).catch(function(e) {});
                    delete cfg.appParams;
                    Ext.apply(app.m_Config, cfg);
                    app.createAndLoadMainPanel();
                }
                else app.loadExtraData();
            }
        }, interval);
    },

    updateTracks: function(tracks) { this.forEachView(function(view) { view.updateTracks(tracks); }); },

    // callback object with success and failure defined
    reset: function(callback, reload) {
        this.m_Config = new SeqView.Config(this);
        if (reload) {
            if (this.m_MarkersInfo) {
                this.m_MarkersInfo.reset();
            }
            this.forEachView( function(view) { view.remove(); } )
            this.m_Views = [];
            this.m_InitialLoading = true;
        }
        var succ_fn = callback && callback.success ? callback.success : function() {};
        var fail_fn = callback && callback.failure ? callback.failure : function() {};
        this.m_Config.load({
            defaultcfg: true,
            resettracks: true,
            callback: {
                success: function(res) { if (reload) { this.infoLoaded(res); } succ_fn(); },
                failure: function(data, text, res) { this.infoMissed(data, text, res); fail_fn(); },
                scope:   this
            }
        });
    },

    // repackage DisplayOptions and TrackConfig to the shape ready for consumption
    // by api functions getDisplayOptions and getTracks
    postProcessConfig: function() {
        // Recode m_Config.DisplayOptions from array to map from track key to track object
        if (!('DisplayOptions' in this.m_Config)) {
            return;
        }
        var tracks = this.m_Config.DisplayOptions.tracks;
        var new_tracks = {};
        for (var i = 0; i < tracks.length; i++) {
            var key = tracks[i].key;
            if ('subkey' in tracks[i]) {
                key += '.' + tracks[i].subkey;
            }
            new_tracks[key] =  tracks[i];
        }
        this.m_Config.DisplayOptions.tracks = new_tracks;
    },

    getTrackObjects: function(all) {
        var tracks = [];
        Ext.each(this.m_Config.TrackConfig, function(track) {
            if (track.shown || all) {
                tracks.push(track);
            }
        });
        tracks.sort(function(t1,t2) {
            return t1.order - t2.order;
        });
        return tracks;
    },

    getTracks: function() {
        return this.getActTracks();
    },

    getActTracks: function() {
         return SeqView.TM.tracksArrayToString(this.getActiveTracks(), true);
    },

    getConfiguration: function() {
        return this.m_Config;
    },

    // Returns an array of track info objects intended for easy constructing
    // 'tracks' URL parameter for subsequence embedding of Sequence Viewer
    // For complete track string either call getTrackStringFromInfo, or embed
    // it's code into your own app.
    getTracksInfo: function() {
        tracks = [];
        for (var i = 0; i < this.m_Config.TrackConfig.length; i++) {
            var track = this.m_Config.TrackConfig[i];
            tracks.push({
                name:  track.name,
                order: track.order,
                shown: track.shown,
                descr: SeqView.TM.trackToString(track)
            });
        }
        return tracks;
    },

    // Returns tracks parameter for Sequence Viewer, use it to construct parameter
    // string for embedding. Does not have parameter name, 'tracks', only value -
    // series of brace enclosed strings describing tracks.
    getTrackStringFromInfo: function(tracks_array) {
        var str = '';
        for (var i = 0; i < tracks_array.length; i++) {
            var track = tracks_array[i];
            str += '[' + track.descr;
            str += ',order:' + track.order;
            str += ',shown:' + (track.shown ? 'true' : 'false');
            str += ']';
        }
        return str;
    },

//////////////////////////////////////////////////////////////////////////
// clearParams:

    clearParams: function(keepViews) {

        this.GI = null;
        this.m_From = null;
        this.m_To = null;
        this.m_ViewSearch = null;
        this.m_ViewSlimModes = [];
        // New markers
        this.m_Markers = [];
        if (!keepViews) {
            this.m_Views = [];
            this.m_Panel = null;
        }
        this.m_MarkersInfo = null;
        this.m_TextView = null;

        this.m_ViewMarkers = false;
        this.m_Key = null;
        this.m_NAA = null;
        this.m_BamPath = null;
        this.m_SRZ = null;
        this.m_hiddenOptions = '';
        this.m_userData = [];
        this.m_userTrackNames = [];
        this.m_FindCompartments = true;
        this.m_showBlastTracks = '';
        this.m_DisplaySettings = '';
        this.m_NeedAlignmentView = false;
        this.m_OnlyAlignmentView = false;
        this.m_DepthLimit = null;
        this.m_AppContext = 'SV-standalone';

        this.m_Origin = 0;
        this.m_Flip = false;
        this.m_ViewLabels = null;
        this.m_Embedded = false;
        this.m_SnpFilter = null;
        this.m_ViewSeq = null;

        this.m_ViewOptions = {};
        this.m_ViewRanges = [];
        this.m_ViewColors = [];
        this.m_ViewsSelect = [];
        this.defaultConfig = {};
        delete this.m_ViewTheme;
        delete this.m_ViewContent;
        delete this.m_AppName;
        delete this.m_AssmContext;

        this.m_ViewParams = null;
        this.m_ItemID = null;

        this.m_LocalPrimerTool = false;
        this.m_LocalBlast = false;
        this.m_ShowQuery = false;
        this.m_QueryRange = null;
        this.m_NoDataCookie = false;
        this.m_SelectionTTMap = new Map();

//        delete this.m_highlightsColor;
        this.m_GraphicExtraParams = {};

        this.m_Toolbar = {
            history: true,
            name: true,
            search: true,
            panning: true,
            zoom: true,
            genemode: true,
            slimmode: true,
            thub: true,
            tools: true,
            download: true,
            config: true,
            reload: true,
            help: true
        };
    },

    parseParams: function(appParams, docParams) {
        delete SeqView.requestTrackSets;

        this.m_AllViewParams = (docParams || '').replace(/\#$/,''); 
        if (!appParams) {
            var rel = Ext.get(this.m_DivId);
            if (rel && rel.dom) {
                var aref = rel.first();
                if (aref && aref.dom && aref.dom.tagName == "A") {
                    appParams = aref.dom.href.slice(aref.dom.href.indexOf('?') + 1)
                    aref.remove();
                } else { // for compatibility with the old verions
                    appParams = rel.dom.getAttribute('rel');
                }
            }
        }
        if (appParams != null) {
            var href = this.m_AllViewParams;
            var keep = ''; 
            var pLst = ['parallel_render', 'extra_opts'];
            Ext.each(pLst, function(p){
                var idx = href.indexOf(p);
                if (idx < 0) return;
                keep += '&' + href.slice(idx).split('&')[0]
            });
            // We should not merge URL parameters in for embedded viewers.
            // Better is to parse appParams into parameter set first,
            // but this is reliable enough - equal sign in all other
            // places beyond parameter definition should be escaped
            this.m_AllViewParams = (/embedded=/.test(appParams))
                ? appParams
                : this.m_AllViewParams + '&' + appParams;
            this.m_AllViewParams += keep;
        }
        this.m_AllViewParams = Ext.util.Format.htmlDecode(this.m_AllViewParams);

        this.m_MarkersInfo = new SeqView.MarkersInfo(this);

        var marker_ranges = [];
        var marker_names = [];
        if (this.m_AllViewParams) {
            var p_array = this.m_AllViewParams.split('&');

            delete this.mf_MultiPanel;
            var sKeys= ['key', 'naa', 'tracks', 'tkey', 'tn', 'rkey', 'mn', 'mk',
                        'url', 'url_reload', 'rid', 'data', 'id', 'select', 'srz',
                        'content', 'viewer_context', 'app_context', 'assm_context',
                        'snp_filter', 'bam_path', 'theme', 'search', 'display_settings'];

            for (var i = 0; i < p_array.length; i++) {
                var pair = p_array[i]
                var the_key = pair.split("=", 1)[0].toLowerCase();
                var val = pair.indexOf('=') > 0 ? unescape(pair.substring(pair.indexOf('=') + 1)) : '';
                
                if (sKeys.indexOf(the_key) == -1) val = val.toLowerCase();
                
                switch (the_key) {
                    case 'bam_path': this.m_BamPath = val; break;
                    case 'base_url': this.m_baseURL = val; break;
                    case 'depth_limit': this.m_DepthLimit = val;  break;
                    case 'id': this.GI = val; break;
                    case 'f': case 'from': this.m_From = val; break;
                    case 't': case 'to':  this.m_To = val; break;
                    case 'gl_debug': if (val == '1' || val == 'true') this.m_GraphicExtraParams['gl_debug'] = 'true'; break;
                    case '_gene': case 'gene': this.m_ViewSearch = val; break;
                    // old way to parse markers
                    case 'm': marker_ranges = val.split(','); break; // like: m=249,276
                    case 'mn': marker_names = val.split(','); break; // like: mn=Marker 1,Marker 2
                    // new way, integrated parameter mk=249:250|Marker 1|ff0000,276|Marker 2|00ff00
                    case 'mk': this.m_MarkersInfo.parseMarkersURL(this.m_Markers, val);  break;
                    case 'key': this.m_Key = val; break;// like: db=nucleotide
                    case 'naa': this.m_NAA = val; break; // like: naa=NA000000003,NA000000004, NA000000004:renderer_name
                    case 'search': this.m_ViewSearch = val; break; // like: search=feature_name
                    case 'srz': this.m_SRZ = val; break; // like: srz=SRZ000200
                    case 'tracks': this.defaultConfig.m_TracksFromURL = this.m_TracksFromURL = val; break;
                    case 'nocache': this.defaultConfig.nocache = val; break;
                    case 'vm': this.m_ViewMarkers = true; break; // like: vm=true
                    case 'tn': this.m_userTrackNames.push(decodeURIComponent(val)); break;
                    case 'find_comp': this.m_FindCompartments = val == "true" || val == "1" || val == "on"; break;
                    case 'show_blast_tracks': this.m_showBlastTracks = val; break;
                    case 'display_settings': this.m_DisplaySettings = val; break;
                    case 'data': val = decodeURIComponent(val); // literal data, URI component encoded
                    case 'rkey':
                    case 'rid':  
                    case 'url_reload': // the same like 'url' with check_cs = true 
                    case 'url': this.m_userData.push([the_key, val]); break; // like: url=www.ncbi.nlm.nih.gov/data.txt
                    case 'align': this.m_NeedAlignmentView = true; break;
                    case 'embedded': this.m_Embedded = val == "true" || val == "1" || val; break;
                    case 'labels': this.m_ViewLabels = val; break;
                    case 'multipanel': this.mf_MultiPanel = (val == '1' || val == 'true' || val == '' || val == null); break;
                    case 'onlyalign': this.m_OnlyAlignmentView = true; break;
                    case 'origin': this.m_Origin = parseInt(val); break;// like: origin=10000 ; 0-based
                    case 'snp_filter': this.m_SnpFilter = val; break;// like: Validated|000010000%20-1%20-1%20-1%20_%20_%20_%2043%20_%20_%20
                    case 'seq': this.m_ViewSeq = val.split(':'); break;// like: seq=2000:7000
                    case 'trackevent': this.m_trackEvent = val == "true" || val == "1" || val; break;// fire event after each track image loaded

//                    case 'highlights_color': this.m_highlightsColor = val; break;// default color for highlights
                    case 'color': // Options are global to SeqView app take only 0-item
                    case 'label': 
                    case 'decor': 
                    case 'spacing': this.m_ViewOptions[the_key] = val.split(',')[0]; break;
                    case 'theme': this.m_ViewTheme = val.split(',')[0]; break; // like: theme=Compact,Details
                    case 'content': this.m_ViewContent = val.split(',')[0]; break;

                    case 'v': if(val.length) { this.m_ViewRanges = val.split(','); } break; // like: v=1000:6000,8000:18723
                    case 'c': this.m_ViewColors = val.split(','); break; // like: color=0000FF,FFFF99
                    case 'slim': this.m_ViewSlimModes = val.split(','); break; // like: slim=1,0,... or slim=true,false,...
                    case 'gflip': this.m_Flip = val=='true' || val=='1'; break; // gflip=true
                    case 'flip': case 'strand': var parts = val.split(','); // like: flip=true,false
                        if (parts && parts.length > 0) this.m_Flip = parts[0]=="true" || parts[0]=="1";
                        break;
                    case 'select': this.m_ViewsSelect = val.split(','); break; // like: select=
                    case 'itemid': this.m_ItemID = val; break;
                    case 'noviewheader': this.m_NoViewHeader = val=='true' || val=='1'; break;
                    case 'viewer_context': this.m_ViewerContext=val; break;
                    case 'app_context': this.m_AppContext=val; break;
                    case 'assm_context': this.m_AssmContext=val; break;
                    case 'noguess_assm': this.m_noGuessAssm = val=='true' || val=='1'; break;
                    case 'nopdf': this.m_NoPDF = val=='true' || val=='1'; break;
                    case 'notrdata': this.m_noTrackData = val=='true' || val=='1'; break;
                    case 'extra_opts': this.m_hiddenOptions = val; break;
                    case 'parallel_render': this.m_parallelRender = val; break;
                    case 'iframe': this.m_iFrame = val; break;
                    case 'nodatacookie': this.m_NoDataCookie = val=='true' || val=='1'; break;
                    case 'showquery': this.m_ShowQuery = val=='true' || val=='1'; break;
                    case 'queryrange': var parts = val.split(':');
                        if (parts && parts.length > 1) {
                            this.m_ShowQuery = true;
                            this.m_QueryRange = [ parseInt(parts[0]), parseInt(parts[1])];
                        }
                        break;

                    case 'primer': if (location.pathname.indexOf('staff') != -1) this.m_LocalPrimerTool = val=='true' || val=='1'; break;
                    case 'blast': if (val == 'web') this.m_LocalBlast = 'web'; else this.m_LocalBlast = val=='true' || val=='1'; break;
                    case 'tkey': this.m_TKey = val; break;
                    case 'appname': if (val && val.length > 0) this.m_AppName = val; break;
                    case 'noconfdlg': this.m_NoConfDlg = val=='true' || val=='1'; break;
                    case 'nosavetracks': this.m_noSaveTracks = true; break;
                    case 'toolbar': if (!val) break;
                        var add = val.charAt(0) != '-';
                        var start = add ? 0 : 1;
                        if (add) this.m_Toolbar = {};
                        for (var j = start; j < val.length; j++) {
                            switch (val.charAt(j)) {
                                case 'b': this.m_Toolbar["history"] = add; break;
                                case 'n': this.m_Toolbar["name"] = add; break;
                                case 's': this.m_Toolbar["search"] = add; break;
                                case 'p': this.m_Toolbar["panning"] = add; break;
                                case 'z': this.m_Toolbar["zoom"] = add; break;
                                case 'g': this.m_Toolbar["genemode"] = add; break;
                                case 'm': this.m_Toolbar["slimmode"] = add; break;
                                case 't': this.m_Toolbar["tools"] = add; break;
                                case 'd': this.m_Toolbar["download"] = add; break;
                                case 'c': this.m_Toolbar["config"] = add; break;
                                case 'r': this.m_Toolbar["reload"] = add; break;
                                case 'h': this.m_Toolbar["help"] = add; break;
                            }
                        }
                        break;
                }
                SeqView.TM.renderStat = this.m_hiddenOptions.indexOf('render_stat') >= 0;
            } // for
            this.m_noTrackData = this.m_Toolbar['download'] != true;
            
            if( this.mf_MultiPanel === undefined ){
                this.mf_MultiPanel = !this.m_Embedded;
            }

            if (this.m_Embedded && !this.m_NoDataCookie) {
                var cookie_name = SeqView.Cookies.UserDataCookieName;
                if (this.m_AppName && this.m_AppName.length > 0 && cookie_name === SeqView.Cookies.UserDataCookieNameBase) {
                    cookie_name += '-' + this.m_AppName;
                }

                if( !this.m_Key ){
                    // if the NC key with the user's data is not passed check if it is set in cookies
                    this.m_Key = SeqView.SessionData.get( cookie_name, null );
                } else {
                    // save it into a session cookie
                    SeqView.SessionData.set( cookie_name, this.m_Key );
                }
            }
        } // view
        this.m_MarkersInfo.parseMarkersOldStyle(this.m_Markers, marker_ranges, marker_names);

        this.m_Portal =
            this.m_Embedded === false
            && (this.m_AppName || '').indexOf('ncbientrez') == 0;

        var tmpl = new Ext.Template(
            '<a onClick="SeqView.App.showLinkURLDlg(' + this.m_Idx + ', \'3-0\');" href="#">Link To This View</a> | ',
            '<a onClick="SeqView.App.showFeedbackDlg(' + this.m_Idx + ', \'3-1\');" href="#">Feedback</a>');
        var cfg = [
            {tag: 'div', cls: 'SeqViewerControls hidden_for_print', id: this.m_DivCtrlId, html:tmpl.apply()},
            {tag: 'div', cls: 'SeqViewerJS', id: this.m_DivJSId}];

        if (this.m_Embedded) {
            cfg.shift();
            cfg[0].cls = '';
        }
        if (this.GI || this.m_userData.length){
            if (!Ext.fly(this.m_DivId).query('div[class="SeqViewerControls hidden_for_print"]').length)
                Ext.DomHelper.append(this.m_DivId, cfg);
            if (!this.m_Panel) {
                this.m_Panel = new Ext.Panel({
                        collapsible: (this.m_Embedded === false),
                        renderTo: this.m_DivJSId,
                        html: '<span id="string_ruler_unit"></span>',
                        header: false
                });
                this.m_Panel.tmpPanel = this.m_Panel.add(new Ext.Panel({
                    height: 60, border: false,
                    items:[{xtype:'displayfield', name: 'progress', id:'sv-uplmsg' + this.m_Idx, value: '',
                             style: {color:'grey', "text-align": 'left', "margin-top":'10px', "margin-left":'6px'}},
                           {xtype:'displayfield', name: 'error', id:'sv-uplerr' + this.m_Idx, value: '',
                             style: {color:'red', "text-align": 'left', "margin-top":'10px', "margin-left":'6px'}}]
                }));
            }
            this.m_Panel.updateLayout();
        }
    },

    parseTracksFromURL: function() {
        return null; // this option is disabled 
        var tParams = this.m_TracksFromURL.replace(/\[/g,'').split(']');
        tParams.pop();
        var trx = [],
            sbt = false;

        for (var i = 0; i < tParams.length; i++) {
            var tp = tParams[i],
                stp;
            if (tp.search(/\\:/) != -1) {
                tp = tp.replace(/:/g, '::').replace(/\\::/g, '\\:')
                stp = tp.split('::')
            } else stp = tp.split(':');
            
            var trk = {cfgIdx: i, tParams: tParams[i]},
                kWord = stp.shift();
            while (stp.length > 1) {
                var v = stp.shift();
                delim = v.lastIndexOf(',');
                trk[kWord] = v.slice(0, delim);
                kWord = v.slice(delim + 1);
            }
            trk[kWord] = stp[0];
            if (!trk.id || trk.id.length < 8 || !trk.key) return null;
            trk.category = { name: trk.category || 'Other' };
            trk.subcategory = { name: trk.subcategory || 'Other' };
            trk.display_name = NCBIGBUtils.unescapeName(trk.display_name || 'No name');
            trk.annots = [NCBIGBUtils.unescapeName(trk.annots || '')];

            if (trk.subtracks) {
                trk.subTracks = trk.subtracks.split('|');
                trk.legend = [];
                sbt = true;
            }
            if (trk.key == 'sequence_track')
                trk.check_boxes = [{value: tParams[i].indexOf(',ShowLabel:false') == -1}];
            trx.push(trk);
        }
        if (!sbt) return trx;
        trx.forEach(function(t) {
            t.shown = t.shown == 'true' || t.shown == "1" || t.shown == 1; 
            if (!t.subTracks) return;
            t.subTracks.forEach(function(stID) {
                var st = trx.find(function(trk) { return trk.id == stID || trk.uid == stID; }) || {};
                st.uId = stID;
                t.legend.push({id: st.id || stID, color: st.color, opasity: st.opacity || 100});
            });
        });
        return trx;
    },

    resizeIFrame: function(height){
        if (!this.m_iFrame) return;
        var h = document.documentElement;
        var height = Math.max(height || 0, h.offsetHeight + (Ext.isIE ? 9 : Ext.isFirefox ? 5 : 0));
        window.frames.parent.postMessage({f: this.m_iFrame, h: height}, '*');
    },

//////////////////////////////////////////////////////////////////////////
// openNewWindowPOST(url, params) - open new window by POSTing 'params' to 'url'
//     url - URL to use for POST
//     params - object with parameters to pass in the POST request
    openNewWindowPOST: function(url, params) {
        var new_win = window.open('about:blank', '_blank');
        new_win.focus();
        var new_body = new_win.document['body'];
        var form = new_win.document.createElement('form');
        form.method = 'POST';
        for (var name in params) {
            if (!params.hasOwnProperty(name)) continue;
            var el = new_win.document.createElement('input');
            el.type = 'hidden';
            el.name = name;
            el.value = params[name];
            form.appendChild(el);
        }
        new_body.appendChild(form);
        form.action = url;
        form.submit();
    },

    blast: function(view, range) {
        var params = {};
        params["SV_SHOW"] = "TRUE";
        if (this.m_Key) {
            params["SV_KEY"] = this.m_Key;
        }
        params["PAGE_TYPE"] = "BlastSearch";
        params["SHOW_DEFAULTS"] = "on";
        if (this.m_ViewParams['acc_type']=='protein') {
            params["PROGRAM"] = "blastp";
            params["BLAST_PROGRAMS"] = "blastp";
        } else {
            params["PAGE"] = "Nucleotides";
            params["PROGRAM"] = "blastn";
            params["MEGABLAST"] = "on";
            params["BLAST_PROGRAMS"] = "megaBlast";
        }
        if (range) {
            params["QUERY_FROM"] = (range[0]+1);
            params["QUERY_TO"]   = (range[1]+1);
        }

        // two next parameters switch BLAST results page to new mode [SV-4937]
        params["ADV_VIEW"] ="on";
        params["SV_SHOW"] = "false";

        var prod = SeqView.webNCBI.indexOf("www.ncbi.nlm.nih.gov") != -1;
        var prefix = this.m_LocalBlast === 'web' ? "https://web.ncbi.nlm.nih.gov/blast/" : (
                                           this.m_LocalBlast ? "" : (prod ? "https://blast.ncbi.nlm.nih.gov/" :
                                           SeqView.webNCBI + "blast/"));
        var url = prefix + "Blast.cgi";
        // check local id, fill out query with sequence in this case
        if (this.GI.indexOf("lcl|") != 0) {
            params["QUERY"] = this.GI;
            this.openNewWindowPOST(url, params);
        } else {
            var ranges = "0-" + (this.m_SeqLength-1);
            var app = this;
            var callback = function(key) {
                var rq = new XMLHttpRequest();
                rq.open('get', app.m_CGIs.NetCache + key);
                rq.onload = function(e) {
                     // The define line is coming with coordinate range, which Blast does not understand.
                    // So, we need to remove it from, e.g. ">lcl|Query_1:1-2534" to read ">lcl|Query_1"
                    params["QUERY"] = rq.response.replace(/>(lcl\|.*):\d+-\d+/, function(m, p1) { return ">" + p1; });
                    params["QUERY_BELIEVE_DEFLINE"]="true";
                    app.openNewWindowPOST(url, params);
                }
                rq.send();
            }
            
            this.requestFile(view, '?id=' + this.GI + '&ranges=' + ranges + '&format=fasta', '', callback);
        }
    },

//////////////////////////////////////////////////////////////////////////
// primerBlast:
    primerBlast: function(whole, range_set) {
        var seqinfo = this.m_ViewParams;
        if (seqinfo["acc_type"] !== "DNA") return;
        if (!whole  &&
            (!range_set || range_set.length < 1 || range_set.length > 2)) return;
        var url = this.m_LocalPrimerTool
            ? "primer-blast/index.cgi?"
            : "/tools/primer-blast/index.cgi?"
        ;
        url += "ORGANISM=" + encodeURI(seqinfo["organism"]);
        url += "&INPUT_SEQUENCE=" + encodeURI(seqinfo["id"]);
        var total_range;
        if (whole)
            total_range = [0, this.m_SeqLength-1];
        else
            total_range = range_set[0];
        if (!whole  &&  range_set.length === 2) {
            var range0, range1;
            if (range_set[0][1] < range_set[1][0]) {
                var range0 = range_set[0];
                var range1 = range_set[1];
            } else {
                var range0 = range_set[1];
                var range1 = range_set[0];
            }
            url += "&PRIMER5_END=" + (range0[1]+1);
            url += "&PRIMER3_START=" + (range1[0]+1);
            total_range = [ range0[0], range1[1] ];
        }
        url += "&PRIMER5_START=" + (total_range[0]+1);
        url += "&PRIMER3_END=" + (total_range[1]+1);
        url += "&PRIMER_PRODUCT_MAX=" + (total_range[1] - total_range[0]);
        url += "&SHOW_SVIEWER=true";
        if (this.m_Key) {
            url += "&SVIEWER_DATA_KEY=" + this.m_Key;
        }
        window.open(url);
    },

    createAndLoadMainPanel: function() {
        this.m_BrowWidth = Ext.getBody().getViewSize().width;
        if (!this.m_Views.length) Ext.fly(this.m_DivId).on({'contextmenu': this.onContextMenu, scope: this});
        this.loadAccession();
    },

    onContextMenu: function(e) {
        if (e.getTarget().id.indexOf('combobox')) e.stopEvent();
    },

    requestGI: function() {
        Ext.get('app-frontpage').setStyle('display', 'block');
        var button = Ext.get('mainshow-btn');
        button.on('click', function() { this.loadAccessionDlg();  },this);
    },

    loadAccessionDlg: function() {
        Ext.MessageBox.prompt('Load Accession', 'Please enter accession or GI:', function(btn, text) {
            if (btn!='ok'  || text.length==0) return;
            this.redirectWithGi(text);
        }, this, false, this.GI);
    },

    redirectWithGi: function(text_gi) {
        var url = this.m_CGIs.Config;
        this.AjaxRequest({url:url  + '?id=' + (text_gi) + '&notrackconfig=true', context: this,
            success:function(data, text, res) {
                var from_cgi = SeqView.decode(data);
                if (from_cgi.job_status) {
                    var st = from_cgi.job_status;
                    if (st == 'submitted' || st == 'running' || st == 'pending') {
                        Ext.defer(SeqView.App.simpleAjaxRequest, 2000, this, {url: url + '?job_key=' + from_cgi.job_id, context: this});
                        return;
                     }
                }
                if (from_cgi.error) {
                    Ext.MessageBox.show({title: 'Sequence Viewer', msg:from_cgi.error,
                                         buttons: Ext.MessageBox.OK, icon:Ext.MessageBox.ERROR});
                } else if(from_cgi.success === false) {
                    Ext.MessageBox.show({title: 'Sequence Viewer', msg:from_cgi.msg,
                                         buttons: Ext.MessageBox.OK, icon:Ext.MessageBox.ERROR});
                } else {
                    var href = document.location.href.replace(/#/, '');
                    var pars = href.split("?");
                    var args = pars[1] ? pars[1].split('&') : [];
                    var view_url = pars[0] + '?id=' + text_gi;
                    Ext.each(args, function(v) {
                        if (v.indexOf('id=') != 0) view_url += '&' + v;
                    });
                    document.location.href = view_url;
                }
            }});
    },

    showTrackHubDlg: function() {
        SeqView.TM.showTrackHubDialog(this);
    },

    showTracksConfigDlg: function(area) {
        this.pingClick(area);
        SeqView.TM.showConfigDialog(this, 0);
    },

    updateConfig: function(newconfig, user_data) {
        this.m_Config.TrackConfig = newconfig.TrackConfig;
        Ext.apply(this.m_Config.Options, newconfig.Options);
        this.reloadAllViews();
        delete this.m_currentTrackSetId;
        this.fireEvent('configuration_changed', this);
    },


//////////////////////////////////////////////////////////////////////////
// showTracks: show/hide tracks matching the query. NB: This function does NOT
//             refresh the view, do it yourself after the call if it returned
//             true.
//    query - object to match with the track to operate on, currently 2 keys
//            supported - "category" and "name". Conditions are joned by AND
//            so if you turn on track with {category: "Variation", name: "SNP" }
//            only tracks with matching category AND names become active.
//            Default (empty object) means ALL tracks.
//    show  - boolean, whether to turn tracks on or off, default off
//    save  - should we save new track configuration, default no
//  returns whether track configuration changed
    showTracks: function(query, show) {
        var active_tracks = [];
        var changed = false;
        Ext.each(this.m_Config.TrackConfig, function(track) {
            var shown = track.shown;
            // concatenate conditions with AND
            var select = false;
            for (var p in query) {
                var pp = {}; 
                for (pp in query[p]) {
                    if (pp === '0') {
                        select |= track[p] === query[p];
                        break;
                    }
                    select |= track[p][pp] === query[p][pp];
                }
             }
            if (select) track.shown = show;
            if (track.shown) active_tracks.push(track);
            if (shown != track.shown) changed = true;
        });
            
        active_tracks.sort(function(t1,t2) { return t1.order - t2.order; });
        return changed;
   },

    getActiveTracks: function(addHidden) {
        var active_tracks = [];
        var subtracks = [];
        var trackConfig = this.m_Config.TrackConfig;
        Ext.each(trackConfig, function(track) {
            if (!track.shown) return true;
            active_tracks.push(track);
            Ext.each(track.subTracks, function(uId, idx) {
                var legend = track.legend[idx];
                if (addHidden || uId.substr(-7) !== '_hidden') subtracks.push(legend);
                Ext.each(trackConfig, function(ct, cfgIdx) {
                    if (ct.id == legend.id) {
                        var opacity = legend.opacity || Math.round(legend.color.split(' ')[3]/2.55);
                        if (!ct.hasOwnProperty('hidden_settings'))
                            ct.hidden_settings = [{name: 'opacity', value: opacity}];
                        delete legend.opacity;
                        Ext.each(ct.hidden_settings, function() { legend[this.name] = this.value; });
                        if (!legend.hasOwnProperty('opacity')) {
                           ct.hidden_settings.push({name: 'opacity', value: opacity});
                           legend.opacity = opacity;
                        }
                        legend.idx = cfgIdx;
                        return false;
                    }
                });
            });
        });
        active_tracks.sort(function(t1,t2) { return t1.order - t2.order; });

        Ext.each(subtracks, function() {
            var trk = trackConfig[this.idx];
            if (trk && !trk.shown) active_tracks.push(trk);
        });
        return active_tracks;
    },

    loadAccession: function() {
        if (!this.GI) return;
        if (!this.m_noTrackData) this.m_noTrackData = this.GI.search(/lcl\||Query_/) == 0;
        if (this.m_MarkersInfo) this.m_MarkersInfo.reset();

//        this.forEachView(function(view) { view.remove(); })
//        this.m_Views = [];
        this.m_InitialLoading = true;

        var cfg = {callback: {success: this.infoLoaded, failure: this.infoMissed, scope: this}};
        if (!this.m_Config.trackInfo && this.m_TKey) {
            var params = {what:'tracks', key: this.m_TKey};
            this.AjaxRequest({
                url: this.m_CGIs.UUD_dlmeta, data: params, context:this,
                success: function(data, text, res){
                    if (data.success === true)
                        this.defaultConfig.m_TracksFromURL = this.m_TracksFromURL = data.tracks;
                    delete this.m_TKey;
                    this.m_Config.load(cfg);
                },
                error: function(data, text, res){
                    delete this.m_TKey;
                    this.m_Config.load(cfg);
                }
            });
        } else {
            this.m_Config.load(cfg);
        }
    },

    adjustViewOptions: function() {
        Ext.each(['label','color','decor','spacing'], function(par) {
            var pidx = this.m_ViewOptions[par] || -1;
            if (pidx >= 0 && pidx < this.m_Config.Options.controls[par].length)
                this.m_Config.Options['curr_' + par] = this.m_Config.Options.controls[par][pidx];
        }, this);
    },

    infoMissed: function(data, text, res) {
        var frontpage = Ext.get('app-frontpage');
        if (frontpage) frontpage.setStyle('display', 'block');
        var button = Ext.get('mainshow-btn');
        if (button) button.on('click', function() { this.loadAccessionDlg(); }, this);

        this.m_InitialLoading = false;
        if (data === undefined) return;
        var dt = SeqView.decode(data || {});
        this.m_Panel.add({html: 'An internal error has occurred that prevents Sequence Viewer from displaying.<br>'
            + 'Technical details (seqconfig error): '
            + (dt.msg || dt.statusText || dt.error_message || text)});
    },

    infoLoaded: function(seqconfig) {
        var frontpage = Ext.get('app-frontpage');
        if (frontpage) frontpage.remove();
        this.m_Config.checkWarnings();
        this.m_ViewParams = this.m_Config.SeqInfo;
        this.m_ViewParams.units = this.m_ViewParams.acc_type == 'protein' ? 'aa' : 'nt';
        this.m_AssmContext = this.m_ViewParams.assm_context;
        this.m_SeqLength = this.m_ViewParams['length'];
        if (!this.m_Config.Options) {
            Ext.applyIf(this.m_ViewOptions, {color: 'Color', decor: 'Default', label: 'Default', spacing: 'Normal', content: 'Details'});
            this.m_Config.Options = {};
            for(var param in this.m_ViewOptions) {
                this.m_Config.Options['curr_' + param] = this.m_ViewOptions[param];
            }
            this.m_Config.Options.controls = this.m_Config.Options.controls || {
                color: ["Color", "Greyscale"],
                decor: ["Default", "Arrows", "CircleAnchor", "Fancy", "SquareAnchor"],
                label: ["Default", "Top label", "Side label", "No label"],
                spacing: ["Normal", "Compact", "Oversize", "Overview"]
            };
        }

        this.postProcessConfig();

        this.fireEvent('configuration_loaded', this);

        if (this.m_Embedded === false)
            document.title = this.m_ViewParams['id'] + ': ' + this.m_ViewParams['title'];

        if (Ext.isFirefox && document.documentElement.style.overflowY == 'scroll'
            && (this.m_Embedded === false || this.m_AppName.search(/demo/i) > 0) ) document.body.style.overflow = 'hidden';

        var title_place = Ext.get(this.m_DivTitle);
        if (title_place) title_place.update(this.m_ViewParams['title']);

        var title_id = Ext.get(this.m_DivTitleID);
        if (title_id) title_id.update(this.m_ViewParams['id_full']);

        this.m_Panel.setTitle(this.m_ViewParams['id']);
        var views = this.m_ViewParams['views'];
        var suggested_ranges = [];
        for (var i = 0; i != views.length; i++) {
            var sp = views[i].split(':');
            // To prevent showing only alignment for PrimerBlast we make an exception here
            if (sp[0]=='multialign' && this.defaultConfig.m_ViewContent != 'PrimerBlast') {
                this.m_NeedAlignmentView = true;
            }
            if (sp[0]=='graphical'  &&  sp[1]) {
                suggested_ranges = sp[1].split('-');
            }
            views[i] = sp[0];
        }

        this.adjustViewOptions();
        if (this.m_Embedded === false || this.m_Embedded == 'panorama' || this.m_Embedded == 'full') {
            this.createPanorama();
        }

        var from = 0,
            to = this.m_SeqLength - 1;
        if (this.m_From && (this.m_From=="begin"  ||  this.m_From=="begining")) this.m_From = null;
        if (this.m_To && this.m_To=="end") this.m_To = null;
        if (!this.m_From && !this.m_To) {
            if (suggested_ranges[0]) {
                from =  NCBIGBUtils.stringToNum(suggested_ranges[0]);
                to  =   NCBIGBUtils.stringToNum(suggested_ranges[1]);
//                m_From = (from + 1) + ''; // convert to string
//                m_To   = (to + 1) + ''; // convert to string
            }
        } else {
            from = this.m_From ? NCBIGBUtils.stringToNum(this.m_From) - 1 : 0;
            to   = this.m_To   ? NCBIGBUtils.stringToNum(this.m_To) - 1   : this.m_SeqLength - 1;
        }

        if (this.m_QueryRange) {
            var overhang = Math.round(0.15 * (this.m_QueryRange[1] - this.m_QueryRange[0]));
            this.m_ViewRanges = [ Math.max(0, this.m_QueryRange[0]-overhang ) + ":" +
                                  Math.min(this.m_SeqLength-1, this.m_QueryRange[1]+overhang)];
        }

        // fix for old Portal format
        if (this.m_ViewRanges.length==1 && this.m_ViewRanges[0].indexOf('begin') != -1){
            this.m_ViewRanges = []; // ignore cases such as begin..end
        }
        if (this.m_Embedded == 'panorama') return;
        if (this.m_ViewRanges.length == 0) {
            this.m_ViewRanges[0] = (from + 1) + ':' + (to + 1);
        }

        if ((this.m_OnlyAlignmentView === true && this.m_Embedded === false) || (this.m_Embedded && this.m_NeedAlignmentView === true)) {
            var alview = new SeqView.Alignment(this);
            var range = this.m_ViewRanges[0].split(/:|\.\.|-/);
            var zs = range[1] == 'zs';
            range = this.decodeRange(range);
            alview.m_FromSeq = range[0];
            alview.m_LenSeq  = this.m_SeqLength;
            if (zs) {
                alview.m_LenSeq -=  Math.floor(415 * SeqView.MinBpp); // adjust to align viewer borders
            }
            this.registerView(alview);
        } else {
            for (var i = 0; i < this.m_ViewRanges.length; ++i) {
                var view = this.m_Views[this.m_Panorama ? i + 1 : i] || new SeqView.Graphic(this);
                var range = this.m_ViewRanges[i].split(/:|\.\.|-/);
                var zs = range[1] == 'zs';

                range = this.decodeRange(range);
                view.m_FromSeq = range[0];
                view.m_LenSeq  = range[1] - range[0] + 1;

                view.m_UrlFrom = range[0] ? range[0] +1 : null;
                view.m_UrlTo   = range[1] ? range[1] +1 : null;
                view.m_slimMode = false;
                var vsm = this.m_ViewSlimModes;
                if (vsm.length > i) view.m_slimMode = (vsm[i] == "true" || vsm[i] == "1");
                if (this.m_ViewColors.length > i) view.m_Color = this.m_ViewColors[i];
                view.setFlipNoReload(this.m_Flip, true);

                if (this.m_ViewsSelect.length > i  && this.m_ViewsSelect[i] != 'null') view.m_SelectedSig = this.m_ViewsSelect[i];

                if (!view.m_View) this.registerView(view);
                else view.createPanel();

                if (i == 0 && this.m_NeedAlignmentView && this.m_Embedded === false) {
                    var alview = new SeqView.Alignment(this);
                    alview.m_FromSeq = range[0];
                    alview.m_LenSeq  = range[1] - range[0] + 1;
                    if (zs) {
                        alview.m_LenSeq -=  Math.floor(415 * SeqView.MinBpp); // adjust to align viewer borders
                    }
                    this.registerView(alview);
                }
            }
        }

        if (this.m_ViewSeq && this.m_ViewSeq.length > 0) {
            var from = NCBIGBUtils.stringToNum(this.m_ViewSeq[0])-1;
            var to = NCBIGBUtils.stringToNum(this.m_ViewSeq[1])-1;
            this.createTextView([from,to,this.m_Flip]);
        }
    },
/*to be removed*/
    getToolTipForTrack: function(view, tname) {
        //for tracks extracting first part, delimited by ";"
        if(tname.match(";")){
            var arr = tname.split(";");
            var gijunk = arr[1];//arr[1] is used for getting track name
            var id = gijunk.split("|");
            tname = id[1];
        }
        return this.m_Config.getTrackByName(tname).help;
    },
/*to be removed*/
    notifyViewLoaded: function(view) {
        this.fireEvent((view.isPanorama() ? 'panorama' : 'graphical') + '_image_loaded', view);
        if (this.m_InitialLoading) {
            var loading_finished = true;
            this.forEachView(function(v) { loading_finished = loading_finished && !v.isLoading(); });
            if (loading_finished) {
                this.m_InitialLoading = false;

                 this.loadSearchPatternData();
                 if (!this.hasMarkers()) {
                    for (var m = 0; m != this.m_Markers.length; m++) {
                        this.addMarkerByData(this.m_Markers[m]);
                    }
                } // create markers from URL
                if (this.m_ViewMarkers && !view.isPanorama()) { // show marker info dialog?
                    this.m_ViewMarkers = false;
                    this.showMarkersDlg(view);
                }
                if (this.m_ViewSearch && this.m_ViewSearch.length > 0) {
                    var search_term = this.m_ViewSearch;
                    this.forEachView(function(v) {
                        var gotoBox = v.m_View.down('#gotoBox');
                        if (gotoBox) gotoBox.setValue(search_term);
                    });
                    this.gotoAndSearch(this.m_ViewSearch);
                }


            }
        }
    },

    moveTo: function(from, len){
        this.forEachView(function(v) {
            v.moveTo(from, len, {from_ui: false});
        });
    },

    addView: function(cfg) {
        var view  = this.m_Panel.add(cfg);
        this.m_Panel.updateLayout();
        return view;
    },

    updateLayout: function() {
        if (this.m_Panel)
            this.m_Panel.updateLayout();
    },

    viewIsClosing: function(view) {
    },

    countActiveSeqPanels: function() {
        var count = 0;
        for (var i = 0; i < this.m_Views.length; i++) {
            if (this.m_Views[i] && this.m_Views[i].m_Type == 'graphical')
                count++;
        }
        return count;
    },

    removeView: function(view) {
        if (this.m_Panorama && view.m_Locator &&  view.m_Locator == this.m_Panorama.m_Locator) {
            this.m_Panorama.m_Locator = this.m_Views[1].m_Locator;
        }
        for (var i = 0; i < this.m_Views.length; i++) {
            var v = this.m_Views[i];
            if (v && v.m_Idx == view.m_Idx) {
                this.m_Views[i] = null;
            }
        }
    },

    createPanorama: function() {
        if (!this.m_Panorama) {
            this.m_Panorama = new SeqView.Panorama(this);
            this.m_Views.unshift(this.m_Panorama);
        }
        this.m_Panorama.updateView();
    },

    loadPanoramaImage: function() {
        if (this.m_Panel && this.m_Panorama) {
            this.m_Panorama.updateView();
        }
    },

    registerView: function(view) {
        this.m_Views.push(view);
        view.createPanel();
        if (NCBIGBUtils.clearAria) NCBIGBUtils.clearAria(Ext.get(this.m_DivId));
    },

    createTextView: function(view) {
        // when the TextView dialog (view) closes it will set this.m_TextView to null
        if (!this.m_TextView)
            this.m_TextView = new SeqView.TextView(this);
        this.m_TextView.openDlg(view);
    },

    getPanoramaHeight: function() {
        return this.m_Panorama? this.m_Panorama.getHeight() : 0;
    },

    getPanoramaWidth: function() {
        return this.m_Panorama? this.m_Panorama.getWidth() : 0;
    },

    updateLocator: function(view) {
        if (view.isPanorama() || view.isLoading()) {
            return;
        }

        var locator = view.m_Locator;
        if (locator && this.m_Panorama) {
            if (!this.m_Panorama.isLoading()) {
                var range = view.toSeq();
                locator.setLeft(this.m_Panorama.toPix(range[0]));
                locator.setWidth(this.m_Panorama.toPix(range[2]));
                locator.setHeight(this.getPanoramaHeight());
            } else {
                Ext.defer(this.updateLocator, 50, this, [view]);
            }
        }
    },

    updateViewTitles: function() {
        this.forEachView(function(view) { view.updateTitle(); });
    },

    findView: function(idx) {
        for (var i = 0; i < this.m_Views.length; i++) {
            var view = this.m_Views[i];
            if (view && view.m_Idx == idx) {
                return view;
            }
        }
        return null;
    },

    forEachView: function(fn, scope) {
        var array = this.m_Views;
        for(var i = 0, len = array.length; i < len; i++){
            var v = array[i];
            if (v) {
                if (fn.call(scope || v, v, i, array) === false) { return (i+1); }
            }
        }
    },

    getFlip: function() {
        return this.m_Flip;
    },

    getOrigin: function() {
        return this.m_Origin;
    },

    setFlip: function(flip, view) {
        this.m_Flip = flip;
        var view_idx = -1;
        if (typeof view != "undefined") {
            view_idx = view.m_Idx;
        }
        this.forEachView(function(view) {
            if (view.isGraphic()  &&  view.m_Idx != view_idx) {
                view.setFlipLocal(flip);
            }
        });
        this.fireEvent('strand_changed', this, view);
    },

    decodeRange: function(range) {
        var from = NCBIGBUtils.stringToNum(range[0])-1;
        var to = 0;
        if (!from) from = 0;
        if (from >= this.m_SeqLength) from = 0;
        if (range[1] == 'zs') {
            var len = Math.floor((this.m_Panel.body.getWidth() - 10) * SeqView.MinBpp);
            to = from+len-1;
        } else {
            to = NCBIGBUtils.stringToNum(range[1])-1;
        }
        if (!to) to = this.m_SeqLength;
        to = Math.min(to, this.m_SeqLength - 1);
        return [from, to];
    },

    // Keep pos or range regex in one place
    splitPosOrRange: function(s) {
        s = s.replace(/[, ]/g, '');
        return s.match(/^([-+]?\d+(?:\.\d+)?[km]?)(?:(-|to|\.\.+|:|\/|_)([-+]?\d+(?:\.\d+)?[km]?))?$/i);
    },

    // s - position or range to handle
    // options - Object with parameters
    //     allow_equal - should we allow both ends of range to be equal, false
    //     ask_user    - if position is ambiguous, can we ask user, false
    //     success     - callback for successfully parsed position
    //         pos_range - array with 1 or two elements representing pos or range
    //                     in sequence coordinates
    //         options   - options passed to handlePos
    //     failure       - callback for reporting errors
    //         error message
    //         options   - options passed to handlePos
    //     scope       - object to use as 'this' for callbacks
    handlePos: function(s, options) {
        var parts = this.splitPosOrRange(s);
        if (parts == null) {
            var hgvs = this.isHGVSExpression(s);
            var snip = s.match(/^([rs]s)([0-9]{3,})(?::.+)?$/);
            var vardb = s.match(/^([en]s)(td|v|sv)([0-9]+)(?::.+)?$/);
            if (hgvs) {
                this.handleHGVS(s, options);
            } else if (snip) {
                this.handleVariations(snip[1]+snip[2], options);
            } else if (vardb) {
                this.handleVariations(vardb[1]+vardb[2] + vardb[3], options);
            } else if (options.failure) {
                options.failure.call(options.scope, "Invalid position/range", options);
            }
            return;
        }
        var pos0 = parts[1];
        var sep = parts[2];
        var pos1 = parts[3];
        if (this.m_Origin && !this.isExplicitPosOrRange(pos0, sep, pos1)) {
            if (options.ask_user) {
                // Construct choices
                var is_range = pos1 && pos1.length > 0;
                var ambiguousItems = [];
                var lim0 = 2, lim1 = 2;
                var num0 = NCBIGBUtils.stringToNum(pos0);
                var num1 = NCBIGBUtils.stringToNum(pos1);
                var fc0 = pos0.charAt(0);
                if (fc0 == '+' || fc0 == '-') lim0 = 1;
                var fc1;
                if (is_range && sep !== '/') fc1 = pos1.charAt(0);
                if (fc1 == null || (fc1 == '+' || fc1 == '-')) lim1 = 1;
                var baseURL = SeqView.base_url;
                for (var i = 0; i < lim0; ++i) {
                    for (var j = 0; j < lim1; ++j) {
                        var item;
                        var flipped = this.getFlip() ? "_flip" : '';
                        if (is_range) {
                            var n0 = (1-i*2)*num0;
                            var n1 = (1-j*2)*num1;
                            if (n1 < n0) {
                                var t = n0; n0 = n1; n1 = t;
                            }
                            var imageSel = 'mp';
                            if (n0 >= 0 && n1 >= 0) {
                                imageSel = "pp";
                            } else if (n0 < 0 && n1 < 0) {
                                imageSel = "mm";
                            } else if (Math.abs(n0) < Math.abs(n1)) {
                                imageSel = "pm";
                            } 
                            if (n0 > 0) n0 = '+' + n0;
                            if (n1 > 0) n1 = '+' + n1;
                            var fullPos = n0 + ' .. ' + n1;
                            var image = '<img aria-label="icon" src="' + baseURL + 'images/choice_range_' + imageSel + flipped + '.png"/>';
                            item = {boxLabel: image + '&nbsp;' + fullPos, name: 'pos', inputValue: fullPos, checked: (i == 0 && j == 0)};
                        } else {
                            var n0 = (1-i*2)*num0;
                            var imageSel = 'm';
                            if (n0 > 0) {
                                n0 = '+' + n0;
                                imageSel = 'p';
                            }
                            var image = '<img aria-label="icon" src="' + baseURL + 'images/choice_pos_' + imageSel + flipped + '.png"/>';
                            item = {boxLabel: image + '&nbsp;' + n0, name: 'pos', inputValue: n0, checked: (i == 0 && j == 0)};
                        }
                        ambiguousItems.push(item);
                    }
                }
                var title = "Ambiguous " + (is_range ? "range" : "position");
                var subtitle = "Choose the correct " + (is_range ? "range" : "position");
                var dlg = new Ext.Window({
                    constrain: true,
                    app: this,
                    modal: true,
                    title: title,
                    width: 400,
                    layout: 'fit',
                    items: [{
                        xtype: 'form',
                        bodyPadding: '5px',
                        items: [
                            {xtype: 'fieldset', title: subtitle, layout: 'form',
                                items: [{ xtype: 'radiogroup', columns: 1, items: ambiguousItems }]},
                            {xtype: 'displayfield',
                                value: 'To avoid ambiguity, the range coordinate values must use explicit + and - signs'
                                + ' when the origin has been reset.'
                        }]
                    }],
                    buttons: [
                        {text: 'OK', scope: this,
                            handler: function() {
                                var values = dlg.items.items[0].getForm().getValues();
                                this.handlePos.call(this, values['pos'], options);
                                dlg.close();
                            }},
                        {text: 'Cancel', handler: function() { dlg.close(); }}
                    ]
                });
//                if (document.body.clientHeight > dlg.height) dlg.renderTo = this.m_DivId;
                dlg.show();
            } else {
                if (options.failure) {
                    options.failure.call(options.scope, "Explicit position/range required", options);
                }
            }
        } else {
            var beg_pos = this.convertRelativePosition(pos0);
            if (pos1 && pos1.length > 0) {
                // range is specified
                var end_pos;
                if (sep == "/") {
                    var pad = NCBIGBUtils.stringToNum(pos1);
                    end_pos = beg_pos + pad;
                    beg_pos -= pad;
                } else {
                    end_pos = this.convertRelativePosition(pos1);
                    if (beg_pos > end_pos) {
                        var t = beg_pos;
                        beg_pos = end_pos;
                        end_pos = t;
                    }
                }
                if (isNaN(beg_pos) || isNaN(end_pos) || beg_pos < 0 || end_pos < 0 ||
                    (options.allow_equal ? beg_pos > end_pos : beg_pos >= end_pos)  ||
                    beg_pos >= this.m_SeqLength || end_pos >= this.m_SeqLength)
                {
                    if (options.failure) {
                        var from_pos = this.posToLocal(0);
                        var to_pos = this.posToLocal(this.m_SeqLength - 1);
                        if (to_pos < from_pos) {
                            var temp = to_pos;
                            to_pos = from_pos;
                            from_pos = temp;
                        }
                        var msg = 'Invalid range: Sequence positions should be from ' + from_pos + ' to ' + to_pos;
                        options.failure.call(options.scope, msg, options);
                    }
                } else {
                    if (options.success) {
                        options.success.call(options.scope, [beg_pos, end_pos], options);
                    }
                }
            } else {
                // single position
                if (isNaN(beg_pos) || beg_pos < 0 || beg_pos >= this.m_SeqLength) {
                    if (options.failure) {
                        var from_pos = this.posToLocal(0);
                        var to_pos = this.posToLocal(this.m_SeqLength - 1);
                        if (to_pos < from_pos) {
                            var temp = to_pos;
                            to_pos = from_pos;
                            from_pos = temp;
                        }
                        var msg = 'Invalid position: Sequence position should be from ' + from_pos + ' to ' + to_pos;
                        options.failure.call(options.scope, msg, options);
                    }
                } else {
                    if (options.success) {
                        options.success.call(options.scope, [beg_pos], options);
                    }
                }
            }
        }
    },

    // term - HGVS position to parse
    // options - Object with parameters
    //     allow_equal - should we allow both ends of range to be equal, false
    //     success     - callback for successfully parsed position
    //         pos_range - array with 1 or two elements representing pos or range
    //                     in sequence coordinates
    //         options   - options passed to handlePos
    //     failure       - callback for reporting errors
    //         error message
    //         options   - options passed to handlePos
    //     scope       - object to use as 'this' for callbacks
    handleHGVS: function(term, options) {
        var view = options ? options.view : null;
        var params = {term: term, id: this.GI, type: "hgvs"};
        if (this.m_Key && this.m_Key.length > 0)
            params.key = this.m_Key;
        var tracks = '';
        this.getActiveTracks().forEach(function(t){
            if (t.key != 'alignment_track') tracks += SeqView.TM.tracksArrayToString([t], true);
        });
        if (tracks)
            params.tracks = tracks;

        return [this.AjaxRequest({
            url: this.m_CGIs.FeatSearch,
            data: params,
            context: this,
            success: function(data, text, res) {
                var from_cgi = SeqView.decode(data);
                if (from_cgi.total) {
                    var feat = from_cgi.features[0];
                    var pos = feat.from;
                    var pos_range = [pos];
                    if (feat.to != pos) {
                        pos_range.push(feat.to);
                    }
                    var label = feat.label;
                    if (options && options.success) {
                        options.success.call(this, pos_range, options, label);
                    } else {
                        this.setPositionalMarker(pos_range, label);
                        if (view) {
                            view.gotoPosRange([pos], true, {from_ui: true} );
                        } else {
                            this.forEachView(function(v) {
                                if (v.isGraphic())
                                    v.gotoPosRange([pos], true, {from_ui: true});
                            });
                        }
                    }
                } else {
                    var msg = from_cgi.errmsg || (term + " not found on the sequence");
                    if (options && options.failure) {
                        options.failure.call(this, msg, options);
                    } else {
                        Ext.MessageBox.alert('Search ' + (from_cgi.errmsg ? 'Error' : 'results'), msg);
                    }
                }
            },
            error: function(data, text, res) {
                var msg = "Server error in parsing/searching request " + term;
                if (options && options.failure) {
                    options.failure.call(this, msg, options);
                }
            }
        })];
    },

    handleVariations: function(term, options) {
        var extSearch = false, found = 0, rqNum = 1,
            pos_range, params,
            trx = this.m_Config.TrackConfig.filter(function(t) { return t.category.name == 'Variation'; });
            app = this,
            view = options ? options.view : null;

        var completeTask = function() {
            if (--rqNum > 0) return;
            if (!found && !extSearch) { //searching in not displayed tracks
                rqNum = 0;
                extSearch = true;
                trx.forEach(function(t, idx) {
                    if (t.shown) return;
                    params.track = '[' + SeqView.TM.trackToString(t) + ']';
                    params.pos = params.from;
                    rqNum++;
                    app.AjaxRequest({url: app.m_CGIs.Graphic, data: params, error: completeTask,
                        success: function(res) { processResponse(res, idx); }});
                });
                return;                
            }
            if (rqNum < 0) return;
            app.unmaskGraphViews(view);
            if (found && extSearch) app.fireEvent('configuration_changed', app);
            if (view) {
                view.gotoPosRange(pos_range, true, {from_ui: true} );
            } else {
                app.forEachView(function(v) {
                   if (v.isGraphic()) v.gotoPosRange(pos_range, true, {from_ui: true});
                });
            }
        }
        var processResponse = function(res, idx) {
            if (res.error_message) {
                app.unmaskGraphViews(view);
                Ext.MessageBox.alert('Search Error', res.error_message);
                return;
            }
            var jst = res.job_status;
            if (jst == 'submitted' || jst == 'running' || jst == 'pending') {
                Ext.defer(SeqView.App.simpleAjaxRequest, 1000 * (rqNum % 5), app,
                    [{url: app.m_CGIs.Graphic, data: {job_key: res.job_id}, error: completeTask, success: function(res) { processResponse(res, idx);}}]);
                return;
            }          
            if (res.new_pos == pos_range[0]) {
                trx[idx].shown = true;
                found++;
                if (!extSearch && rqNum > 0) rqNum = 1;
            }
            completeTask();
        }

        var params = {term: term, id: this.GI, type: "snp"};
        if (this.m_Key) params.key = this.m_Key;
        app.maskGraphViews(view, 'Searching...');
        return [this.AjaxRequest({
            url: this.m_CGIs.FeatSearch,
            data: params,
            context: this,
            success: function(data, text, res) {
                var from_cgi = SeqView.decode(data);
                if (from_cgi.total) {
                    var snp_found = false;
                    for (var i = 0; i < from_cgi.features.length; i++) {
                        var feat = from_cgi.features[i];
                        if (feat.seqId != this.GI) continue;
                        snp_found = true;
                        var pos = feat.from;
                        pos_range = [pos];
                        if (feat.to != pos) pos_range.push(feat.to);
                        var label = feat.label;
                        if (options && options.success) {
                            options.success.call(this, pos_range, options, label);
                        } else {
                            this.setPositionalMarker(pos_range, label, view);
                            params = {id: app.GI, navi: 1, len: 10, width: 100, from: pos + (!pos ? 1 : -1), dir: (!pos) ? 'prev' : 'next'};
                            // searching in displayed tracks
                            trx.forEach(function(t, idx) {
                                if (!t.shown) return;
                                params.track = '[' + SeqView.TM.trackToString(t) + ']';
                                params.pos = params.from;
                                rqNum++;
                                app.AjaxRequest({url: app.m_CGIs.Graphic, data: params, error: completeTask,
                                    success: function(res) { processResponse(res, idx); }});
                            });
                            completeTask();
                        }
                        break;
                    }
                    // it's a valid snp but it's not annotated on this sequence
                    // so we suggest exploring dbSNP page
                    if (!snp_found) {
                        app.unmaskGraphViews(view);
                        var snip = term.match(/^([rs]s)([0-9]{3,})$/);
                        var snip_url = SeqView.webNCBI + 'snp/' + snip[2];
                        Ext.Msg.alert("SNP search", 'Valid SNP is not found on this sequence, please see <a target="_blank" href="'+snip_url+'">dbSNP page for '+term+'</a>');
                    }
                } else {
                    app.unmaskGraphViews(view);
                    var msg = term + " not found on the sequence";
                    if (from_cgi.errmsg) console.log(msg, from_cgi.errmsg);
                    if (options && options.failure) {
                        options.failure.call(this, msg, options);
                    }
                }
            },
            error: function(data, text, res) {
                app.unmaskGraphViews(view);
                var msg = "Server error in parsing/searching request " + term;
                if (options && options.failure) {
                    options.failure.call(this, msg, options);
                }
            }
        })];
    },

// parsePosOrRange:
// returns an array of one or two elements, depending on input, which
//    represent a valid position or range on sequence. If there is no valid
//    position/range can be parsed, returns undefined. All input is in local
//    coordinate system - 1-based with flip and origin taken into account.
//    Output is 0-based, in sequence coordinates.

    parsePosOrRange: function(s, allow_equal) {
        s = s.replace(/[, ]/g, '');
        var parts = this.splitPosOrRange(s);
        if (parts == null) {
            return;
        }
        var part0 = parts[1];
        var sep = parts[2];
        var part1 = parts[3];
        var beg_pos = this.convertRelativePosition(part0);
        if (part1 && part1.length > 0) {
            // range is specified
            var end_pos;
            if (sep == "/") {
                var pad = NCBIGBUtils.stringToNum(part1);
                end_pos = beg_pos + pad;
                beg_pos -= pad;
            } else {
                end_pos = this.convertRelativePosition(part1);
                if (this.getFlip()) {
                    var t = beg_pos;
                    beg_pos = end_pos;
                    end_pos = t;
                }
            }
            if (isNaN(beg_pos) || isNaN(end_pos) || beg_pos < 0 || end_pos < 0 ||
                (allow_equal ? beg_pos > end_pos : beg_pos >= end_pos)  ||
                beg_pos >= this.m_SeqLength || end_pos >= this.m_SeqLength)
            {
                return;
            }
            return [beg_pos, end_pos]
        } else {
            // single position
            if (isNaN(beg_pos) || beg_pos < 0 || beg_pos >= this.m_SeqLength)
                return;
            return [beg_pos];
        }
    },


    isPosOrRange: function(s) {
        return this.splitPosOrRange(s) != null;
    },

//////////////////////////////////////////////////////////////////////////
// isExplicitPosOrRange:
//     pos0, sep, pos1 - parts of a range split by pattern
//     returns true if numerals have a sign. e.g. -12:+24, or +35..-20
    isExplicitPosOrRange: function(pos0, sep, pos1) {
        var fc0 = pos0.charAt(0);
        var fc1;
        if (pos1 && pos1.length > 0 && sep !== '/') fc1 = pos1.charAt(0);
        return !((fc0 != '+' && fc0 != '-') || (fc1 && fc1 != '+' && fc1 != '-'));
    },

    convertRelativePosition: function(pos_str) {
        var pos = NCBIGBUtils.stringToNum(pos_str);
        if (!isNaN(pos))
            return this.posToGlobal(pos);
        return pos;
    },

//////////////////////////////////////////////////////////////////////////
// posToLocal: convert backend 0-based global coordinate to
//     local, taking into account strand and origin
    posToLocal: function(pos, flip) {
        pos -= this.m_Origin;
        if (pos >= 0) pos += 1;
        return pos;
    },

//////////////////////////////////////////////////////////////////////////
// posToLocalDisplay: convert backend 0-based global coordinate to
//     local, taking into account origin, returns number with sign if
//     origin is set, so for positive numbers it will return a string
    posToLocalDisplay: function(pos, flip) {
        pos -= this.m_Origin;
        if (pos >= 0) {
            pos += 1;
            if (this.m_Origin) pos = '+' + pos;
        }
        return pos;
    },

//////////////////////////////////////////////////////////////////////////
// posToGlobal: convert 1-based coordinate relative to origin and strand
//     to global 0-based backend coordinate
    posToGlobal: function(pos) {
        if (pos > 0) pos -= 1;
        return pos + this.m_Origin;
    },

    newMarkerDlg: function(view, x_pos) {
        if (!this.m_MarkersInfo) {
            this.m_MarkersInfo = new SeqView.MarkersInfo(this);
        }
        this.m_MarkersInfo.MarkerDlg(view || this, x_pos);
    },

    showMarkersDlg: function( view ){
        if( !this.m_MarkersInfo ){
            this.m_MarkersInfo = new SeqView.MarkersInfo( this );
        }
        this.m_MarkersInfo.showDlg( view );
    },

    updateMarkersSize: function(view) {
        if (this.m_MarkersInfo) {
            this.m_MarkersInfo.updateMarkersSize(view);
        }
    },

    hasMarkers: function() {
        return this.m_MarkersInfo && this.m_MarkersInfo.hasMarkers();
    },

    getMarkersInfo: function() {
        return this.m_MarkersInfo;
    },

    updateMarkersPos: function(view) {
        if (this.m_MarkersInfo) {
            this.m_MarkersInfo.updateMarkersPos(view);
        }
    },

    forEachMarker: function(fn, scope) {
        if (this.m_MarkersInfo) {
            this.m_MarkersInfo.forEachMarker(fn,scope);
        }
    },

    findMarker: function(m_name) {
        if (this.m_MarkersInfo) {
            return this.m_MarkersInfo.findMarker(m_name);
        }
        return null;
    },

    addMarker: function(marker_data) {
        if (!this.m_MarkersInfo) {
            this.m_MarkersInfo = new SeqView.MarkersInfo(this);

        }
        this.m_MarkersInfo.addMarker(marker_data);
    },

    addMarkerByData: function(marker_data) {
        if (!this.m_MarkersInfo) {
            this.m_MarkersInfo = new SeqView.MarkersInfo(this);
        }
        this.m_MarkersInfo.addMarkerByData(marker_data);
    },

    scrollMarkers: function(view, delta) {
        if (this.m_MarkersInfo) {
            this.m_MarkersInfo.scrollMarkers(view,delta);
        }
    },


    reloadAllViews: function(options) {
        this.forEachView(function(view) { view.refresh(options); });
    },

    showOriginDlg: function(pos) {
        Ext.MessageBox.prompt('Set Origin', 'Please enter new origin:', function(btn, text) {
            if (btn!='ok'  || text.length==0) return;
            var position = 0;
            var bad_pos = !NCBIGBUtils.isNumeric(text);
            if (!bad_pos) {
                position = NCBIGBUtils.stringToNum(text);
                if (position > this.m_SeqLength) {
                    bad_pos = true;
                }
            }
            if (bad_pos) {
                Ext.MessageBox.show({title: 'Set Origin',msg: 'Invalid sequence origin.',
                                     buttons: Ext.MessageBox.OK,icon:Ext.MessageBox.ERROR});
                return;
            }
            position -= 1;
            if (this.m_Origin == position)
                return;
            this.m_Origin = position;
            if (this.m_Origin < 0) this.m_Origin = 0;
            // reload all views
            this.reloadAllViews();
            this.fireEvent('origin_changed', this);
        }, this, false, pos ? pos+1 : this.m_Origin+1);
    },

    setOrigin: function(view, x_pos) {
        var pix = x_pos - view.m_ScrollPix;
        var seq_pos = view.pix2Seq(pix);
        this.showOriginDlg(seq_pos);
    },

    clearOrigin: function() {
        Ext.MessageBox.confirm('Confirm', 'Reset Sequence Origin?', function(btn) {
            if (btn!='yes' || this.m_Origin == 0) return;
            this.m_Origin = 0;
            // reload all views
            this.reloadAllViews();
            this.fireEvent('origin_changed', this);
        }, this);
    },

    loadExtraData: function() {
        this.load_params = this.load_params || {};
        if (!this.m_userData.length) {
            if (this.GI) this.createAndLoadMainPanel();
            else this.infoMissed();
            return;
        }

        if (typeof this.m_AssmContext == 'undefined') {
            this.m_AssmContext = '';
            if (this.GI && this.GI.indexOf("lcl|") !=0) {
                var options = { context: this,
                    url: this.m_CGIs.Config + '?notrackconfig=true&id=' + this.GI,
                    success: function(data, txt, rq) {
                        if (data.job_status) {
                            var st = data.job_status;
                            if (st == 'submitted' || st == 'running' || st == 'pending') {
                                options.url = this.m_CGIs.Config + '?job_key=' + data.job_id;
                                Ext.defer(SeqView.App.simpleAjaxRequest, 1000, this, [options]);
                                return;
                            }
                       }
                       if (data.SeqInfo) this.m_AssmContext = data.SeqInfo.assm_context;
                       this.loadExtraData();
                    },
                    error: this.loadExtraData
                };
                this.AjaxRequest(options);
                return;
            }
        }

        var config = {assm_acc: this.m_AssmContext};//, data_action:'uploading', ret_prj_key: 'true' };
        if (this.GI) config['accession'] = this.GI;
        if (!(config.prj_key = this.m_Key || this.load_params.key)) delete config.prj_key;
        var data_item = this.m_userData.shift();
        if (data_item[0] == 'rkey') {// Fetch long RID with filters from NetCache key
            var params = {data_action: 'downloading', format: 'rids', key: data_item[1], fmt: "text/plain"};
            var processResponse = function(data) {
                    if (data.status == 200) this.m_userData.unshift(['rid', data.responseText]);
                    else  this.m_userTrackNames.shift(); // Keep name array in sync
                    this.loadExtraData();
            }
            this.AjaxRequest({url: this.m_CGIs.NetCache, data: params, context: this,
                success: processResponse,
                error: processResponse});
        } else {
            config.track_name = this.m_userTrackNames.shift();
            config.display_settings = this.m_DisplaySettings;
            switch (data_item[0]) {
                case 'rid': config.blast = {rid: data_item[1], link_related_hits: this.m_FindCompartments, show_blast_tracks: this.m_showBlastTracks};
                    this.load_params.rid_loaded = true; // Mark that we loaded RID
                    break;
                case 'data':
                    config.data = data_item[1];
                    config.file_format = 'asn text';
                    break;
                case 'url_reload': config.check_cs = true; // no break;
                case 'url': config.dataURL = data_item[1]; break;
            }
            var msg = 'Uploading your data';
            if (this.m_AssmContext) {
                msg += " on assembly " + this.m_AssmContext;
            }
            msg += ', please wait...'; 
            this.showMessage(msg);
            this.m_Panel.items.items[0].getEl().mask('Uploading \"' + data_item[0] + '\"');
            var app = this;
            app.consError = console.error;
            console.error = function() {}
            try {
                var uploaderUUD = SeqView.App.getFileUploader(config);
                var promise = uploaderUUD.getPromise();
                promise.fail(function() {
                    console.error = app.consError || console.error; delete app.consError;
                    app.m_Panel.items.items[0].getEl().unmask();
                    var errMsg = this.getErrors();
                    app.showMessage(errMsg, true);
                    var fbText = 'SViewer initial parameters: ' + app.m_AllViewParams + '\nError message: ' + this.getErrors();
                    Ext.Msg.show({title: 'User data loading error',
                        msg: errMsg,
                        buttons: Ext.Msg.YESNOCANCEL,
                        buttonText: {yes: 'Continue', no: 'Feedback'},
                        icon: Ext.MessageBox.WARNING,
                        fn: function(btn) {
                            switch (btn) {
                                case 'no': app.showFeedbackDlg('2-Feedback', fbText, document.location.href); // "break" is missed to continue data loading
                                case 'yes': app.loadExtraData(); break;
                            }
                        }});
                });
                promise.done(function(tlist, dkey) {
                    console.error = app.consError || console.error; delete app.consError;
                    if (app.load_params.rid_loaded && !app.GI) app.initRID(this.getMetadata());
                    app.addKey(dkey);
                    app.addUUDTracks(this.getTracks());
                    var cookie_name = SeqView.Cookies.UserDataCookieName;
                    if (cookie_name === SeqView.Cookies.UserDataCookieNameBase) {
                        cookie_name += (app.m_AppName ? ('-' + app.m_AppName) : '');
                    }
                    SeqView.SessionData.set(cookie_name, app.m_Key);//save the key into a session cookie
                    var msg = 'Data uploaded';
                    if (app.m_AssmContext) 
                        msg += " on assembly " + app.m_AssmContext;
                    app.showMessage(msg);
                    app.loadExtraData();
                });
                var currTask = 'Uploading your data'
                promise.progress(function(progress) {
                    var task = progress.current_task;
                    if (task == "" || task == 'pending' || task == currTask) return;
                    app.showMessage(currTask = task);
                });
                uploaderUUD.upload();
            } catch(e) { app.showMessage('Unable to upload data: ' + e.message, true); } 
        }
    },

    addUUDTracks: function(trx) {
        this.uudTracks = this.uudTracks || [];
        trx.forEach(function(trk) {
            trk.id = trk.id || trk.GetTMSId();
            if (!this.uudTracks.some(function (t) { return t.id == trk.id; })) this.uudTracks.push(trk);
        }, this);
    },

    showMessage: function(msg, errorFlag) {
        var msg_field = Ext.getCmp((errorFlag ? 'sv-uplerr' : 'sv-uplmsg') + this.m_Idx);
        if (msg_field) {
            msg_field.show();
            msg_field.update(msg);
        }
    },

    maskGraphViews: function(view, msg, msgCls) {
        var mask = function(v) { if (v.isGraphic()) Ext.get(v.m_DivId).mask(msg, msgCls); }
        if (view) mask(view); else this.forEachView(mask);
    },
    unmaskGraphViews: function(view) {
        var unmask = function(v) { if (v.isGraphic()) Ext.get(v.m_DivId).unmask(); }
        if (view) unmask(view); else this.forEachView(unmask);
    },

    addKey: function(key) {
        if (typeof key == 'undefined') return;
        if (!this.m_Key || this.m_Key.length == 0) this.m_Key = key;
        else
            if (this.m_Key.search(key) == -1) this.m_Key += ',' + key;
    },

    initRID: function(from_cgi) {
        if (from_cgi['blast_query'] === undefined) return; // Can't do anything
        var blast_query_list = from_cgi['blast_query'];
        var blast_query = {};
        for (var i = 0; i < blast_query_list.length; i++) {
            var kv = blast_query_list[i];
            blast_query[kv.key] = kv.value;
        }
        this.GI = blast_query['id'];
        if (!this.m_QueryRange) {
            var beg = blast_query['beg'];
            var end = blast_query['end'];
            var querybeg = blast_query['querybeg'];
            var queryend = blast_query['queryend'];
            if (querybeg !== undefined  &&  queryend !== undefined) {
                // Convert to numbers
                querybeg = querybeg - 0;
                queryend = queryend - 0;

                // Adjust viewrange, add special marker to show query if requested
                beg = beg !== undefined ? Math.min(beg, querybeg) : querybeg;
                end = end !== undefined ? Math.min(end, queryend) : queryend;
                if (this.m_ShowQuery) {
                     this.m_Markers.push([ [querybeg, queryend], "Query", SeqView.MarkerFlags.SystemLock | SeqView.MarkerFlags.Hollow, "red"]);
                }
            }
            if (beg !== undefined  &&  end !== undefined) {
                var overhang = 0.15 * (end - beg);
                beg = Math.max(0, beg - overhang);
                end += overhang;
                beg = Math.round(beg);
                end = Math.round(end);
                this.m_ViewRanges = [ beg + ':' + end ];
            }
        }
    },

    onSyncAlignView: function(view) {
        for (var i = 0; i < this.m_Views.length; ++i) {
            aview = this.m_Views[i];
            if (aview && aview.isAlignment() ) {
                aview.loadImage(view.m_VisFromSeq, view.m_VisLenSeq);
                break;
            }
        }
    },

    getGraphicViews: function() {
        var views = []
        this.forEachView(function(view) {
            if (view.isGraphic())
                views.push(view);
        });
        return views;
    },

    AjaxRequestExt: function(cfg){
       Ext.applyIf (cfg.params, {appname: this.m_AppName});
       return SeqView.App.simpleAjaxRequestExt(cfg);
    },

    AjaxRequest: function(cfg){
/*
        if( !this.m_AjaxTrans ) this.m_AjaxTrans = [];

        cfg.callback = cfg.callback || function(arg) {};
        cfg.params = cfg.params || {};
*/
        Ext.applyIf (cfg.data, {appname: this.m_AppName});

/*
        var trans = this.m_AjaxTrans;

        var prefunction = function(opts, succsess, res) {
            if( opts.params.transId ){
                trans.remove( opts.params.transId );
                delete opts.params.transId;
            }
        }
        cfg.callback = Ext.Function.createInterceptor(cfg.callback, prefunction);
*/
       return SeqView.App.simpleAjaxRequest(cfg);

    },

//    In a better world IE 11 would support default parameters
//    GridAjaxRequest: function(cfg, req_delay=2000) {
    GridAjaxRequest: function(cfg, req_delay) {
        if (req_delay === undefined) req_delay = 2000;
        function checkJobStatus(data, text, res) {
            var context = cfg.context ? cfg.context : this;
            var from_cgi = SeqView.decode(data);
            if (from_cgi.job_status) {
                if (from_cgi.job_status == 'failed') {
                    if (cfg.error) {
                        cfg.error.call(context, null, from_cgi.error_message);
                    }
                } else if(from_cgi.job_status == 'canceled') {
                    if (cfg.error) {
                        cfg.error.call(context, null, 'Job canceled');
                    }
                } else {
                    var new_cfg = {};
                    Object.assign(new_cfg, cfg);
                    Object.assign(new_cfg, {
                        data: {job_key:from_cgi.job_id},
                        success: checkJobStatus});
                    Ext.defer(this.AjaxRequest, req_delay, this, [new_cfg]);
                }
            } else {
                if (from_cgi.success === false) {
                    if (cfg.error) cfg.error.call(context, null, from_cgi.msg);
                } else {
                    if (cfg.success) cfg.success.call(context, data, text, res);
                }
            }
        }
        var new_cfg = {};
        Object.assign(new_cfg, cfg);
        Object.assign(new_cfg, {success: checkJobStatus});
        this.AjaxRequest(new_cfg);
    },

    getCustomToolTipTools: function(selection) {
        if (this.m_CustomSelectionHandler)
            return this.m_CustomSelectionHandler.getToolTipTools(selection);
        return null;
    },
    setCustomSelectionHandler: function(handler) {
        this.m_CustomSelectionHandler = handler;
    },

    addCustomFeatureFlags: function(cfg) {
        if (this.m_CustomSelectionHandler) {
            this.m_CustomSelectionHandler.addCustomFeatureFlags(cfg);
        }
    },

    setTooltipPreprocessor: function(callback) {
        if (typeof callback === 'function') this.m_preprocessorTT = callback;
        else delete this.m_preprocessorTT;
    },

    watchdogStart: function(url, job_id, params, callback) {
        this.watchdogStop();
        this.m_WatchUrl = url;
        this.m_WatchJobId = job_id;
        this.m_WatchJobParams = params;
        this.m_WatchJobCallback = callback || Ext.emptyFn;
        this.m_WatchTimeoutId = Ext.defer(this.watchdogReport, 90000, this, ['No response']);
    },

    // submit error info to applog
    submitError: function(reason, url, reqparams, jobid) {
        this.AjaxRequest({
            url: this.m_CGIs.Watchdog,
            data: {
                reason: reason,
                requrl: url || this.m_WatchUrl,
                reqparams: reqparams || this.m_WatchJobParams,
                jobid: this.m_WatchJobId
            }
        });        
    },

    watchdogStop: function() {
        if (this.m_WatchTimeoutId ){
            clearTimeout(this.m_WatchTimeoutId);
            this.m_WatchTimeoutId = null;
        }
    },

    watchdogReport: function(reason) {
        this.watchdogStop();
        this.m_WatchJobCallback(); //Ext.defer(this.m_WatchJobCallback(), 30000);
        this.submitError(reason);
    }
});

Ext.apply(SeqView.App, {

    findAppByDivId: function(div_id) {
        for (var i = 0; i < SeqView.m_Apps.length; ++i) {
            if (SeqView.m_Apps[i].m_DivId == div_id) return SeqView.m_Apps[i];
        }
        return null;
    },

    findAppByIndex: function(app_idx) {
        for (var i = 0; i < SeqView.m_Apps.length; ++i) {
            if (SeqView.m_Apps[i].m_Idx == app_idx) return SeqView.m_Apps[i];
        }
        return null;
    },

    getApps: function() { return SeqView.m_Apps; },

    getFileUploader(config) {
        TMS.seqid = config.hint_acc;
        return new UUD.FileUploader(config);
    },

    simpleAjaxRequestExt: function(cfg) {
        cfg.withCredentials = true;
        Ext.applyIf(cfg, {method: 'POST'});
        cfg.cors = true;
//        Ext.applyIf(cfg, {dataType: SeqView.jsonType});
        try {
            return Ext.Ajax.request(cfg);
        } catch (e) { return e;};
    },

    simpleAjaxRequest: function(cfg) {
        cfg.xhrFields = {withCredentials: true};
        Ext.applyIf(cfg, {type: 'POST'});
        cfg.crossDomain = true;
        Ext.applyIf(cfg, {dataType: SeqView.jsonType});
        jQuery.support.cors = true;
        try {
            return jQuery.ajax(cfg);
        } catch (e) { return e;};
    },

    showLinkURLDlg: function(app_idx, logarea) {
        var app = this.findAppByIndex(app_idx);
        if (app) { app.showLinkURLDlg(logarea); }
    },

    showFeedbackDlg: function(app_idx, logarea) {
        var app = this.findAppByIndex(app_idx);
        if (app) app.showFeedbackDlg(logarea);
    }
});

// Portal integration function. It is called then Portal needs to know current
// sequence viewer position.
SeqView.PortalSeqGraphicsInfo = function() {
    var params = "from=&to=&itemid=&strand=";
    if(SeqView.App) {
        var apps = SeqView.App.getApps();
        if (apps && apps.length == 1) {
            var app = apps[0];
            var views = app.getGraphicViews();
            if (views && views.length == 1) {
                var view = views[0];
                if (view.m_VisFromSeq >= 0 && view.m_VisLenSeq && view.m_VisLenSeq > 0) {
                    params = 'from='+(view.m_VisFromSeq+1)+'&to='+(view.m_VisFromSeq+view.m_VisLenSeq);
                } else {
                    params = 'from=&to=';
                }
                params += '&itemid=';
                if (app.m_ItemID) {
                    params += app.m_ItemID;
                }
                if (view.canFlip()) {
                    if (view.getFlip()) {
                        params += '&strand=true';
                    } else {
                        params += '&strand=false';
                    }
                } else {
                    params += '&strand=';
                }
            }
        }
    }
    return params;
};

/*  $Id: sviewapp_more.js 47222 2022-11-15 21:09:42Z lotovv $
 * Authors:  Vlad Lebedev, Maxim Didenko, Victor Joukov
 * File Description:
 */

Ext.apply(SeqView.App.prototype, {

    eventHandlers: function(fn, obj) {
       // fn - 'on'|'un'
       obj = obj.scope ? obj : { scope: obj, pointermove: obj.onPointerMove, pointerup: obj.onPointerUp };
       this.m_Panel.el[fn](obj);
    },
    
    addURLParams: function(params) {
        var app = this.m_App;
        params = params || {};
        if (this.m_NAA) params.naa = this.m_NAA;
        if (this.m_BamPath) params.bam_path = this.m_BamPath;
        if (this.m_SRZ) params.srz = this.m_SRZ;
        if (this.m_Key) params.key = this.m_Key;
        if (this.m_Origin != 0 ) params.origin = this.m_Origin;
        if (this.m_DepthLimit) params.depthlimit = this.m_DepthLimit;
        if (this.m_AssmContext) params.assm_context = this.m_AssmContext;
//        if (this.m_highlightsColor) params.highlights_color = this.m_highlightsColor;

        if (this.m_Config) {
            Ext.apply(params, this.m_Config.visualOptionsUrl())
        }
        return params;
    },

    showLinkURLDlg: function (logarea) {
        if (logarea) this.pingClick(logarea);
        this.getLinkToThisPageURL('', null, function (linkUrl, longUrl) {
            this.resizeIFrame(500);
            longUrl = longUrl || linkUrl;
            var longParams = longUrl.substr(longUrl.indexOf('?') + 1) + '&appname=no_appname';
            var templates = {
                iframe: '<iframe id="iframe_@sv_id@" width="' + this.m_Views[0].m_Width + '" src="' + SeqView.base_url
                      + 'embedded_iframe.html?iframe=iframe_@sv_id@&' + longParams + '" onload="'
                      + 'if(!window._SViFrame){_SViFrame=true;window.addEventListener(\'message\','
                      + 'function(e){if(e.origin==\'https://' + document.domain + '\' && !isNaN(e.data.h))'
                      + 'document.getElementById(e.data.f).height=parseInt(e.data.h);});}">\n</iframe>',
                div: '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">\n'
                      + '<html>\n<head>\n<title>Embedded Sequence Viewer with parameters</title>\n'
                      + '<script type="text/javascript" src="https://www.ncbi.nlm.nih.gov/projects/sviewer/js/sviewer.js"></script>'
                      + '</head>\n<body>'
                      + '<div id="sviewer_@sv_id@" class="SeqViewerApp" data-autoload>\n<a href="?embedded=true&' + longParams + '"></a>\n</div>\n</body>\n</html>'};

            var warnMsg = {title: 'Warning(s)', icon: Ext.Msg.WARNING, buttons: Ext.Msg.OK,
                msg: !this.m_Key ? '' : 'This page contains a link to the user loaded data which is kept<br>' + 
                     'inside a temporary storage and will not be available in approximately 1-2 months.<br>'};

           if (longUrl && longUrl != linkUrl)
                warnMsg.msg += 'The tracks configuration string for this link is too long and it has been saved<br>' +
                           'into a temporary storage. It will be kept there for 3 months.';

            var showEmbedCode = function (tname) {
                var uname = linkDlg.down('#tinyURL').getValue().toLowerCase();
                var parts = uname.split('/');
                if (parts.length > 1) uname = parts.pop().trim();
                else uname = 'NCBI';
                linkDlg.down('#embedCode').setValue(templates[tname].replace(/@sv_id@/g, uname));
            }
            var app = this;
            var linkDlg = new Ext.Window({
                app: this,
                layout: 'vbox',
                modal: true,
                title: 'Link To This View',
                width: 700,// autoHeight: true,
                height: 450,
                constrain: true,
                cls: 'SeqViewerApp',
                listeners: {
                    'beforeshow': function (qt) {
                        NCBIGBUtils.makeTinyURL(linkUrl, function (res) {
                            linkDlg.down('#tinyURL').setValue(res);
                            showEmbedCode('iframe');
                        });
                    }
                },
                resizable: true,
                closeAction: 'destroy',
                plain: true,
                bodyStyle: 'padding:0px;',
                border: false,
                labelAlign: 'right',
                items: [
                    { xtype: 'textarea', height: 115, allowBlank: true, width: '100%', flex:1, value: linkUrl },
                    { xtype: 'textfield', itemId: 'tinyURL', hidden:true, allowBlank: true, width: '100%', emptyText: 'Loading...', value: '' },
                    { xtype: 'fieldset',  width: "100%", flex:2,
                      autoHeight: true,
                      title: 'Embedding code',
                      layout: 'vbox', // 'form',
                      items: [
                        { xtype: 'radiogroup', width: '100%',
                            items: [
                              { boxLabel: "Embed&nbsp;code&nbsp;for&nbsp;IFRAME", name:'rb-embed', checked: true },
                              { boxLabel: "Whole&nbsp;page&nbsp;example", name:'rb-embed',
                                  listeners: {
                                      change: function(fld, newval) { showEmbedCode(!newval ? 'iframe' : 'div'); }
                                  }}
                            ]},
                        { xtype: 'textarea', itemId: 'embedCode', allowBlank: true, flex:1, width: '100%', emptyText: 'Loading...', value: '' }
                    ]}
                ],
                buttons: [!warnMsg.msg ? '' : { itemId: 'warn', xtype: 'tbtext', minWidth: 30, height: 19, text: '', qtip: 'Tooltip', tooltip: 'Warning(s)',
                            style: 'background:url(' + SeqView.ExtIconLoading.slice(5, -18) + 'shared/warning.gif) no-repeat 6px 2px;',
                            handler: function () { Ext.Msg.show(warnMsg); }},
                        '->', { text: 'Close', handler: function () { linkDlg.close(); app.resizeIFrame(); } }]
            });
            linkDlg.show();
            if (warnMsg.msg) {
                var w = linkDlg.down('#warn');
                w.tooltip = Ext.create('Ext.tip.ToolTip', {target: w.id, title: 'Warning(s)'});
                w.tooltip.setHtml(warnMsg.msg);
            }
            
        });
    },


    getLinkToThisPageURL: function (mode, tracks, callback) {
        callback = callback || function (link) { window.open(link); };
        var info = this.m_Config.SeqInfo;
        if (mode == 'GDV') {
            return SeqView.webNCBI + '/genome/gdv/browser/genome/?id=' + info.assm_context + '&chr=' + info.id
                + this.m_MarkersInfo.getMarkersURL();
        }
        tracks = tracks || SeqView.TM.tracksArrayToString(this.getActiveTracks(true), true, true, {tmsSearch: true});
        var appParams = this.getAppParams();
        var link = SeqView.base_url + '?id=' + this.GI;
        if ((this.m_Portal || !mode.indexOf('portal')) && mode != 'full' && info.gi) {
            var db = info.acc_type == 'DNA' ? 'nuccore' : 'protein';
            link = SeqView.webNCBI + db + '/' + info.gi + '?report=graph';
            if (mode == 'portalGB') {
                return link.replace(/graph$/, db == 'nuccore' ? 'genbank' : 'genpept'); 
            }
        }
        var longUrl = link + '&tracks=' + tracks + appParams;
        if (longUrl.length < 2048)
           callback.call(this, longUrl); 
        else {
            var params = { data_action: 'uploading', input_form: 'tracks', tracks: tracks };
            this.GridAjaxRequest({ url: this.m_CGIs.UUD_svdl, data: params, context: this,
                success: function (data, text, res) {
                    callback.call(this, data.success ? link + '&tkey=' + data.key + appParams : longUrl, longUrl);
                },
                error: function (data) { callback.call(this, longUrl); }
            });
        }
    },
    
    
    getAppParams: function(noRange) {
        var params = '';
        Ext.each(
            [['m_Origin','origin'],
             ['m_NAA', 'naa'], ['m_SRZ', 'srz'],
             ['m_BamPath', 'bam_path'],
             ['m_DepthLimit', 'depthlimit'],
             ['m_SnpFilter', 'snp_filter'],
             ['m_ViewerContext','viewer_context'],
             ['m_Key', 'key'],
             ['m_AssmContext', 'assm_context']],
            function(pair) { if (this[pair[0]]) params += '&' + pair[1] + '=' + this[pair[0]]; },
            this);
        if (this.m_AppContext && this.m_AppContext != 'SV-standalone') params += '&app_context=' + this.m_AppContext;
        params += this.m_MarkersInfo.getMarkersURL();

        params += (this.m_noGuessAssm ? '&noguess_assm=1' : '');
        // add views
        var view = '&v';
        var color = '&c';
        var slimMode = '&slim';
        var select = '&select';
        this.forEachView(function (v) {
            if (!v.isGraphic()) return;
            if (v.m_UrlFrom && v.m_UrlTo) { // if we have presise positions from URL, use them. Otherwise - calculate from visible range
                view += ',' + v.m_UrlFrom + ':' + v.m_UrlTo;
            } else {
                view += ',' + (v.m_VisFromSeq + 1) + ':' + (v.m_VisFromSeq + v.m_VisLenSeq);
            }
            color += ',' + v.m_Color;
            slimMode += ',' + (0 + v.m_slimMode);
            select += ',' + (v.m_SelectedSig ? v.m_SelectedSig : 'null');
        });

        var options = this.m_Config.Options;
        if (options) {
            Ext.each(['color', 'label', 'decor', 'spacing'], function (attr) {
                var idx = options.controls[attr].indexOf(options['curr_' + attr]);
                if (idx > 0) params += '&' + attr + '=' + idx;
            });
        }
        params += (view.replace(/,/, '=') + color.replace(/,/, '=') + (this.m_Flip ? '&gflip=1' : '')
                 + select.replace(/,/, '=') + slimMode.replace(/,/, '='));
        params += this.m_Toolbar['download'] ? '' : '&toolbar=-d';
        params += this.m_TextView ? this.m_TextView.getURL() : '';

        if (this.m_SearchDlg && this.m_SearchDlg.isVisible()) {
            var view = null;
            this.forEachView(function (v) { if (!view && v.isGraphic()) view = v; });
            if (view) {
                var combo = view.m_View.down('#gotoBox');
                var searchTerm = null;
                if (combo) searchTerm = combo.getValue();
                if (searchTerm) params += '&search=' + searchTerm.replace(/&/g, '%26');
            }
        }
        return params;
    },

    //////////////////////////////////////////////////////////////////////////
    // gotoAndSearch - if term is HGVS expression, SNP, or position - find the
    // position requested. If it's search term - conduct the search.
    // parameters:
    // term - term to search or go to
    // options - Object with following possible fields
    //     view - if set, update this particular view, otherwise all views
    //     success - if set, report results through this callback,
    //               if not for position set it to the view, or all views if
    //               view is not specified. Set positional marker.
    //               For search, show search dialog.
    //         results - object of the following composition
    //             type - "hgvs",
    //             total - number of features found
    //             features - array of objects
    //                 from - first feature coordinate
    //                 to   - second feature coordinate
    //                 gi   - gi feture found on
    //                 label - feature label
    //         options - options Object passed to gotoAndSearch
    //     failure - if set, report errors through this callback
    //               if not, shows Message Boxes
    //         err_msg - error message text
    //         options - options Object passed to gotoAndSearch
    // returns array with deferred ids
    gotoAndSearch: function (term, options) {
        if (!term) {
            return;
        }
        term = Ext.util.Format.trim(term);
        if (!term.length) {
            return;
        }
        if (!options) {
            // Guard against empty options
            options = {};
        }
        // Add view and scope to callback so it works with parseAndGotoPosition
        if (!options.scope) options.scope = this;
        var hgvs = this.isHGVSExpression(term);
        var snip = term.match(/^([rs]s)([0-9]{3,})(?::.+)?$/);
        var vardb = term.match(/^([en]s)(td|v|sv)([0-9]+)(?::.+)?$/);
        if (hgvs) {
            return this.handleHGVS(term, options);
        } else if (snip) {
            return this.handleVariations(snip[1] + snip[2], options);
        } else if (vardb) {
            return this.handleVariations(vardb[1] + vardb[2] + vardb[3], options);
        } else if (this.isPosOrRange(term)) {
            return this.parseAndGotoPosition(term, options);
        } else {
            if (!options.success) {
                this.openSearchDlg(term, options.view);
                return [];
            } else {
                return this.startSearch(term, options);
            }
        }
    },


    isHGVSExpression: function (term) {
        return term.match(/^(?:([^:]+):)?([cgmnrp])\.(.+)$/);
    },

    //////////////////////////////////////////////////////////////////////////
    // setPositionalMarker - set a marker with a special property: initially locked
    // when unlocked it changes its title to default "Marker #"
    setPositionalMarker: function (pos_range, name, view) {
        if (view) view.zoomSeq(pos_range[0]);
        if (!this.getMarkersInfo().findMarkerByName(name)
            && !(pos_range.length == 1 && this.getMarkersInfo().findMarkerByPos(pos_range[0]))
        ) {
            this.getMarkersInfo().addMarker(pos_range, name,
                SeqView.MarkerFlags.PositionTitle | SeqView.MarkerFlags.UserLock);
        }
    },




    //
    // term - position or range to parse and goto
    // options - Object with parameters:
    //     range_only - should we accept strictly range or position is OK, false
    //     view       - goto this position in a specific view, otherwise - globally
    //     success    - do not go to position, report successfully parsed position to this callback
    //         results in 'features' format
    //         options - options object for parseAndGotoPosition
    //     failure     - do not report error to user, use this callback
    //         error message
    //         options - options object for parseAndGotoPosition
    //     scope      - object to use as 'this' for callbacks
    parseAndGotoPosition: function (term, options) {
        term = term.replace(/[, ]/g, '');
        this.handlePos(term, {
            // pass parameters from our options
            range_only: options.range_only,
            view: options.view,
            nested_success: options.success,
            nested_failure: options.failure,
            // handlePos parameters
            ask_user: true,
            success: function (pos_range, options) {
                if (options.range_only && pos_range.length !== 2) {
                    options.failure.call(options.scope, "Range required", options);
                    return;
                }
                if (options.nested_callback) {
                    var pos = pos_range[0];
                    var pos_to = pos;
                    if (pos_range.length > 1) pos_to = pos_range[1];
                    var res = { type: "pos", total: 1,
                        features: [{ from: pos,
                            to: pos_to,
                            gi: options.app.m_ViewParams.gi,
                            label: term
                        }]
                    };
                    options.nested_callback.call(options.scope, res, options);
                } else {
                    if (pos_range.length == 1) {
                        options.app.setPositionalMarker(pos_range, term, options.view);
                    }
                    if (options.view) {
                        options.view.gotoPosRange(pos_range, true, { from_ui: true });
                    } else {
                        options.app.forEachView(function (v) {
                            if (v.isGraphic()) v.gotoPosRange(pos_range, true, { from_ui: true });
                        });
                    }
                }
            },
            failure: function (err_msg, options) {
                if (options.nested_failure) {
                    options.nested_failure.call(options.scope, err_msg, options);
                } else {
                    Ext.MessageBox.alert('Error', err_msg);
                }
            },
            scope: options.scope ? options.scope : this,
            app: this
        });
        return [];
    },

    checkSearchStatus: function (data, request_type, request_term, options) {
        var from_cgi = SeqView.decode(data);
        if (from_cgi.job_status) {
            var st = from_cgi.job_status;
            if (st == 'submitted' || st == 'running' || st == 'pending') {
                Ext.defer(SeqView.App.simpleAjaxRequest, 2000, this, [{
                    url: this.m_CGIs.FeatSearch + '?job_key=' + from_cgi.job_id,
                    context: this,
                    success: function (data, text, res) {
                        this.checkSearchStatus(data, request_type, request_term, options);
                    },
                    error: function (data, text, res) {
                        this.reportSearchFailure(request_type, request_term, options);
                    }
                }]);
                return;
            }
        }
        from_cgi.type = request_type;
        if (options.success)
            options.success.call(this, from_cgi);

    },
    reportSearchFailure: function (request_type, request_term, options) {
        var msg = "Server error in searching request " + request_term;
        if (options.failure)
            options.failure.call(this, msg);
    },

    startSearch: function (term, options) {
        var seq_type = 'nucleotide';
        if (this.m_ViewParams['acc_type'] == 'protein') {
            seq_type = 'protein';
        }
        var data_key = (this.m_App.m_Key && this.m_App.m_Key.length > 0) ? this.m_App.m_Key : '';
        var tracks = '';
        this.m_App.getActiveTracks().forEach(function(t){
           if (t.key != 'alignment_track') tracks += SeqView.TM.tracksArrayToString([t], true);
        });

        var requests = [
            { url: this.m_CGIs.FeatSearch,
                params: { term: term, type: 'feature', id: this.GI, whole: 1, limit: 1000000, key: data_key, tracks: tracks }
            },
            { url: this.m_CGIs.FeatSearch,
                params: { term: term, type: 'component', id: this.GI, whole: 1, limit: 1000000, key: data_key }
            },
            { url: this.m_CGIs.FeatSearch,
                params: { term: term, type: seq_type, id: this.GI, whole: 1, limit: 1000000, key: data_key }
            }
        ];
        var res = [];

        for (var i = 0, l = requests.length; i < l; i++) {
            res.push((function (that, type, term) {
                return that.AjaxRequest({
                    url: requests[i].url,
                    data: requests[i].params,
                    context: that,
                    success: function (data, text, res) {
                        that.checkSearchStatus(data, type, term, options);
                    },
                    error: function (data, text, res) {
                        that.reportSearchFailure(type, term, options);
                    }
                });
            })(this, requests[i].params.type, requests[i].params.term));
        }
        return res;
    },

    /// MRU search pattern list support


    spsRecId: 100,
    searchPatternData: [],

    loadSearchPatternData: function () {
        var app = this;
        this.AjaxRequest({
            url: this.m_CGIs.SearchMru, context: this,
            data: {
                session: SeqView.SessionData.get("ncbi_sid") || '',
                seqid: this.GI,
                oper: 'load'
            },
            success: function (data, text, req) {
                if (data.mrulist) {
                    var mrulist = Ext.urlDecode(data.mrulist).list;
                    if (mrulist) {
                        Ext.each(Ext.util.JSON.decode(mrulist), function(pt) {app.addSearchPattern(pt, 'push');});
                        app.forEachView(function() {
                            var combo = this.m_View.down('#gotoBox');
                            if (combo) combo.setStore(app.searchPatternData);
                        });
                    }
                }
            },
            error: function () { }
        });

    },

    addSearchPattern: function(pattern, func) {
        if (typeof pattern != 'string') return;
        var maxNumOfSearchPatterns = 20;
        var spd = this.searchPatternData;
        var ix = spd.indexOf(pattern);
        if (ix >= 0) spd.splice(ix, 1); 
        spd[func || 'unshift'](pattern);
        spd.splice(maxNumOfSearchPatterns);
    },

    saveSearchPatternData: function () {
        var mrulist = Ext.urlEncode({list:  Ext.util.JSON.encode(this.searchPatternData)});
        this.AjaxRequest({
            url: this.m_CGIs.SearchMru, context: this,
            data: {
                session: SeqView.SessionData.get("ncbi_sid") || '',
                seqid: this.GI,
                oper: 'save',
                mrulist: mrulist
            },
            success: function (data, text, req) { },
            error: function (data, text, req) { }
        });
    },

    showSearchParamsDlg: function (view) {
        var combo = Ext.create('Ext.form.field.ComboBox', {
            listeners: {
                specialkey: function(f,e){
                    if (e.getKey() == e.ENTER){
                       this.searchRequest();
                    }
                }
            },
            typeAhead: true,
            store: this.searchPatternData,
            
           searchRequest: function() { 
               var term = combo.getValue();
               if (!term) return;
               dlg.close();
               view.gotoAndSearch(term);
           }
        });
        var dlg = new Ext.Window({
            title: 'Find on Sequence',
            constrain: true,
            app: this,
            modal: true,
            layout: 'form',
            width: 420,
            padding: 3,
            buttonAlign: 'center',
            buttons: [
                {text: 'OK', handler: combo.searchRequest},
                {text: 'Cancel', handler: function(){ dlg.close(); }}]
        });
        dlg.add({xtype: 'displayfield', bodyStyle: {background: 'inherit'}, border: false,
            value: 'Enter sequence position or range (possible range formats are '
            + '10k-20k, -20--10, -10k:-5, 5 to 515, -1m..1m)<br/><br/>'
            + 'Or enter name of feature, component, HGVS, SNP rs id, '
            + 'or sequence (nucleotide regexp with IUPAC equivalents, PROSITE patterns):<br/>'})

        dlg.add(combo);
        dlg.show();
    },
/*
    showSearchDlg: function () {
        var str = this.m_View.down('#gotoBox').getValue();
        openSearchDlg(str);
    },
*/
    openSearchDlg: function (str, view) {
        var str = str || '*';
        if (this.m_SearchDlg) {
            this.m_SearchDlg.close();
            this.m_SearchDlg = null;
        }
        var fixDelay = function(ajax, res, opt) {
            var from_cgi = res.responseType != 'json' ? Ext.JSON.decode(res.responseText) : res.responseJson;
            opt.callback = opt.callback || opt.callbackBAK;
            if (from_cgi.job_status && from_cgi.job_id) {
                opt.callbackBAK = opt.callback || opt.callbackBAK;
                delete opt.callback;
                opt.params = {job_key: from_cgi.job_id};
                Ext.defer(Ext.Ajax.request, 2000, this, [opt]);
           } 
           return true;
        }
                
        if (!Ext.Ajax.hasListeners.requestcomplete) {
            Ext.Ajax.on('requestcomplete', fixDelay) 
        }
        
        var searchTab = function (title, view, str, tabpanel) {
            this.m_View = view;
            var params = { term: str, id: this.m_View.m_App.GI, from: this.m_View.m_FromSeq,
                to: (this.m_View.m_FromSeq + this.m_View.m_LenSeq - 1),
                limit: 50, whole: 1, appname: this.m_View.m_App.m_AppName
            };
            if (this.m_View.m_App.m_Key && this.m_View.m_App.m_Key.length > 0) {
                params.key = this.m_View.m_App.m_Key;
            }
            switch (title) {
                case "Components":
                    params.type = 'component';
                    break;
                case "Sequence":
                    if (view.m_App.m_ViewParams['acc_type'] == 'protein') {
                        params.type = 'protein';
                    } else {
                        params.type = 'nucleotide';
                    }
                    break;
                case "Features":
                    params.type = 'feature';
                    var tracks = '';
                    view.m_App.getActiveTracks().forEach(function(t){
                        if (t.key != 'alignment_track') tracks += SeqView.TM.tracksArrayToString([t], true);
                    });
                    if (tracks)
                        params.tracks = tracks;
                    break;
            }
            var proxy = {
                type: 'ajax',
                url: view.m_App.m_CGIs.FeatSearch, extraParams: params,
                reader: { rootProperty: 'features', totalProperty: 'total' },
                simpleSortMode: true
            };
            var fields = [
                'label',
                {name: 'from', type: 'int'},
                {name: 'to', type: 'int'},
                {name: 'strand', mapping: function(data)
                    {return data.strand == '-' ? 'Negative' : 'Positive';}}
            ];
            


            var store = Ext.create('Ext.data.Store', {
                model: Ext.define(null, {
                    extend: 'Ext.data.Model',
                    fields: fields,
                    proxy: proxy }),
                pageSize: 50,
                remoteSort: true,
                autoLoad: true,
                job_id: null,
                sorters: [{ property: 'label', direction: 'ASC' }],
                listeners: {
                    load: function (obj, data) {
                        if (this.panel.isDestroyed) return;
                        this.panel.setIconCls(data && data.length > 0 ? 'xsv-search-results' : '');
                }}
            });
            store.proxy.setActionMethods({read: 'POST'});
            var panel = new Ext.grid.GridPanel({
                title: title,
                store: store,
                iconCls: 'xsv-search-loading',
                border: false,
//                loadMask: true,
                stripeRows: true,
                columns:[
                    { text: "Label", flex: true, wi_dth: 110, sortable: true, dataIndex: 'label' },
                    { text: "From", width: 90, sortable: true, dataIndex: 'from', renderer: function (val) { return val + 1; } },
                    { text: "To", width: 90, sortable: true, dataIndex: 'to', renderer: function (val) { return val + 1; } },
                    { text: "Strand", width: 90, sortable: true, dataIndex: 'strand'} 
                ],
                viewConfig: {forceFit:true, deferEmptyText:false, emptyText:'<div align="center">No Search Results To Display</div>'},
                listeners: {
                    'rowdblclick': function (table, rec, tr, rowIndex, e) {
                        var view = this.m_View;
                        if (table.ownerGrid.title == 'Component') {
                            var tracks_changed = view.m_App.showTracks({ key: 'component_track' }, true);
                            tracks_changed |= view.m_App.showTracks({ key: 'scaffold_track' }, true);
                            if (tracks_changed) view.m_App.fireEvent('configuration_changed', view.m_App);
                        }
                        var from = rec.get('from');
                        var to = rec.get('to');
                        var view_len = to - from + 1;

                        view.m_SelectedSig = rec.get('object_id') || null;
                        if (view.m_SelectedSig)
                            view.removeRangeSelection(true);
                        else
                            view.m_SelectedRangeSet =[[from, to]];
                        view.pingClick('2-3-' + table.ownerGrid.title.charAt(0));
                        view.startImageLoading(from - view_len * .15, view_len * 1.3);
                        view.m_App.fireEvent('selection_changed', view); 
                    },
                    scope: this
                },
                bbar: new Ext.PagingToolbar({
                    store: store,
                    displayInfo: true, displayMsg: 'Displaying Search Results {0} - {1} of {2}', emptyMsg: "No Search Results To Display"
                })

            });
            store.panel = panel;
            panel.on('beforedestroy', function () {
                if (this.store.job_id != null) {
                    SeqView.App.simpleAjaxRequest({ url: this.store.url + '?Cancel=true&job_key=' + this.store.job_id });
                }
            });
            return panel;
        };

        // This is a hack. we need to find a first graphic view in order to get the search range
        // ideally, the search function should be called with the view parameter
        this.forEachView(function (v) { if (!view && v.isGraphic()) view = v; });
        if (view) {
            var tabPanel = new Ext.TabPanel({
                activeTab: 0,
                listeners: {
                    'tabchange': function (tabpanel, tab) {
                        tabpanel.down('#reconf_btn').setVisible(tab.title == 'Tracks');
                    }
                },
                buttons: [
                    {text: 'Re-configure', hidden: true, scope: this, itemId: 'reconf_btn',
                        handler: function () {
                            this.pingClick('2-3-T');
                            if (this.m_SearchDlg.getActiveTab().syncModifiedTracks()) {
                                var categories = SeqView.TM.processTracksInfo(this.m_Config.TrackConfig);
                                SeqView.TM.Common.updateSeqViewApp(categories, this);
                            }
                            this.m_SearchDlg.close(); this.m_SearchDlg = null;
                        }
                    },
                    {text: 'Close', scope: this, handler: function () { this.m_SearchDlg.close(); this.m_SearchDlg = null; } }
                ]
            });
            tabPanel.add(new searchTab("Features", view, str, tabPanel));
            tabPanel.add(new searchTab("Components", view, str, tabPanel));
            tabPanel.add(new searchTab("Sequence", view, str, tabPanel));
            tabPanel.add(SeqView.TM.createSearchTracksTab("Tracks", this.m_Config.TrackConfig, str));
            this.m_SearchDlg = new Ext.Window({
                app: this,
                layout: 'fit',
                title: 'Search Results',
                minWidth: 480, width: 650, height: 350,
                constrain: true,
                collapsible: true,
                cls: 'SeqViewerApp',
                closeAction: 'destroy',
                plain: true,
                items: [tabPanel],
                getActiveTab: function ()
                { return tabPanel.activeTab; }
            });
            var app = this;
            this.m_SearchDlg.on('close', function () {
                app.resizeIFrame(); app.m_DialogShown = false;
            });

            app.resizeIFrame(400);
            app.m_DialogShown = true;
            this.m_SearchDlg.show();
        }
    },

    requestFile: function(view, params, fname, callback, nomask) {
        var url = this.m_CGIs.SequenceSave;
        var processResponse = function(from_cgi, t, r) {
            if (from_cgi.job_status) {
                if (from_cgi.job_status == 'failed' || from_cgi.job_status == 'canceled') {
                    Ext.Msg.show({ buttons: Ext.MessageBox.OK, icon: Ext.Msg.ERROR,
                        title: 'Job ' + from_cgi.job_status,
                        msg: from_cgi.error_message });
                    if (!nomask) this.m_App.unmaskGraphViews(view);
                } else {
                    Ext.defer(SeqView.App.simpleAjaxRequest, 500, this, [{
                        url: url + '?job_key=' + from_cgi.job_id, context: this,
                        success: processResponse, error: processResponse }]);
                }
                return;
            }
            if (!nomask) this.m_App.unmaskGraphViews(view);
            if (from_cgi.success === true) {
                if(from_cgi.file_size === 0) {
                    Ext.MessageBox.show({ title: 'Error: ', msg: 'No file created', buttons: Ext.MessageBox.OK, icon: Ext.Msg.ERROR });
                } else {
                    if (callback) {
                        callback(from_cgi.file_url);
                        return;
                    }
                    var form = Ext.DomHelper.append(document.body,
                        { tag : 'form', method : 'post', action : this.m_App.m_CGIs.NetCache + from_cgi.file_url + fname });
                    document.body.appendChild(form);
                    form.submit();
                    document.body.removeChild(form);
                }
            } else
                Ext.MessageBox.show({ title: 'Error: ', msg: from_cgi.msg || t, buttons: Ext.MessageBox.OK, icon: Ext.Msg.ERROR });
        }

        if (this.m_Key) params += '&key=' + this.m_Key;
//        if(params.indexOf('format=fasta') > 0) {
//            url = url.replace('sequence.cgi', 'fasta_sequence.cgi'); 
//        }
        this.AjaxRequest({url: url + params + '&netcache=1', context: view, success: processResponse, error: processResponse});
        if (!nomask) this.maskGraphViews(view, 'Preparing...');
    },

    showFeedbackDlg: function (logarea, fbText, fbURL)
    {
        this.pingClick(logarea);



        var app = this;
        if (typeof fbText != 'string') fbText = null;
        var fbCallback = function (linkUrl)
        {
            let encodedLinkUrl = encodeURIComponent(linkUrl);

            NCBIGBUtils.makeTinyURL(encodedLinkUrl, function (res)
            {
                let encodedTinyUrl = (res!="")?encodeURIComponent(res):encodedLinkUrl;
                window.open("https://nlmenterprise.co1.qualtrics.com/jfe/form/SV_5iL6W3r2A7vY52R?ncbi_app=SV&CurrentPage="+encodedTinyUrl, '_blank');
            });




        }
        if (fbText) fbCallback(fbURL);
        else this.getLinkToThisPageURL('', null, fbCallback);
    },

    findSelectionToolTip: function(key) { return this.m_SelectionTTMap.get(key); },

    saveSelectionToolTip: function(key, value) { this.m_SelectionTTMap.set(key, value); },

    completeCfgLoading: function(cfg) {
        this.forEachView(function(v) { v.updateStatusInfo(); if (!v.isGraphic()) return; v.syncTracks(cfg.TrackConfig); });
    }


});

//////////////////////////////////////////////////////////////////////////
// SeqView.PrinterFriendlyPage


SeqView.PrinterFriendlyPage = function(app, printdlg) {
    this.values = {orientation: 1, sequence: true, title: true}
    var img_width = this.values['orientation']==2 ? 900 : 650;
    var marker_num = 0;
    this.generator=window.open('','name','width='+(img_width+50)+',height=700,toolbar=0,location=0,status=0,menubar=1,resizable=1,scrollbars=1');
    this.generator.document.write('<html><head><title>' + app.GI + '</title>');
    this.generator.document.write('<link rel="stylesheet" href="' + app.m_CGIs.baseURL + 'css/style.css">');
    var callFixSeqTrans = false;
    if (this.values['sequence'] && app.m_TextView && app.m_TextView.m_HideTrans) {
        this.generator.document.write('<style>');
        if (app.m_TextView.m_HideTrans & 1) {
            this.generator.document.write('.xsv-seqtrans_trans {display:none;}');
            this.generator.document.write('.xsv-seqtrans_trans_row {display:none;}');
        } 
        if (app.m_TextView.m_HideTrans & 2) {
            this.generator.document.write('.xsv-seqtrans_prot {display:none;}');
            this.generator.document.write('.xsv-seqtrans_prot_row {display:none;}');
        }
        this.generator.document.write('</style>');
    }
    this.generator.document.write('</head><body onLoad="window.print()">');
    if (this.values['title'])  {
        this.generator.document.write('<p><b>'+app.m_ViewParams['title']+'</b></p>'); // add title
        this.generator.document.write('<p>'+app.m_ViewParams['id_full']+'</p>'); // add full id
    }

    this.ondataloadcalled = false;
    this.running_requests = 0;
    this.overview_data = null;

    var view_data_idx = 0;
    this.views_data = [];

    this.text_data = null;
    if (this.values['sequence'] && app.m_TextView) { // checkbox checked
        this.col_tmp = app.m_TextView.m_SequenceCols;
        app.m_TextView.m_SequenceCols = this.values['orientation'] == 2 ? (Ext.isWindows ? 100: 110) : (Ext.isWindows ? 70 :80);
        if( app.m_TextView.m_Flip ){
            var numLines = app.m_TextView.m_LenSeq / app.m_TextView.m_SequenceCols;
            numLines = Math.ceil( numLines );
            app.m_TextView.m_LenSeq = numLines * app.m_TextView.m_SequenceCols;
        }
        var the_url = app.m_TextView.getSeqTextURL();
        
        this.running_requests++;
        app.AjaxRequest({url: the_url, context: this,
            success:function(data, text, res) {
                this.text_data = SeqView.decode(data);
                this.running_requests--;
                this.onDataLoaded();
            }
       });
    }

    if ( this.ondataloadcalled === false && this.running_requests == 0) {
        // only sequence title is printed
        this.generator.document.write('</body></html>');
        this.generator.document.close();
    }

    this.onDataLoaded = function() {
        this.ondataloadcalled = true;

        if (this.running_requests != 0)
            return;
        if (this.overview_data) {
            // If the img_url begins with ? it contains only parameters for ncfetch, so prepend ncfetch URL
            // This is a way to provide reliable URL resolution for embedding. SV-1760
            var img_url = this.overview_data.img_url;
            if (img_url && img_url.charAt(0) == '?') {
                img_url = app.m_CGIs.NetCache + img_url;
            }
            this.generator.document.write('<img aria-label="icon" style="border:1px solid;" src="' + img_url + '">');
        }
        var print_tpl_str = '<div style="position:absolute;top:2px;width:16px;height:16px;left:{left}px;"><div class="marker_label" style="{font_size}border:1px {color} solid;color:black;">{trimmed_label}</div><div class="marker_line" style="border-right:1px solid {color}; height:{height}px;"></div></div>'
        Ext.each(this.views_data, function(view_data) {
            var view = view_data.view;
            var from_cgi = view_data.data;
            this.generator.document.write('<br><br>'+view.m_View.title + '<br>');
            var html = '<div style="position:relative;">';
            var margin = 0;
            if (app.m_MarkersInfo) {
                app.m_MarkersInfo.forEachMarker( function(the_m) {
                    if ( the_m.seq_pos>view.m_VisFromSeq && the_m.seq_pos<(view.m_VisFromSeq+view.m_VisLenSeq-1)) {
                        var bpp = from_cgi.len / img_width;
                        var marker_left = 0;
                        if (view.getFlip())
                            marker_left = (from_cgi.len - (the_m.seq_pos - from_cgi.from) -1 + 0.5) / bpp - 8;
                        else
                            marker_left = (the_m.seq_pos + 0.5 - from_cgi.from) / bpp - 8;

                        var options = {color:the_m.color, left:marker_left, height:from_cgi.h + 4,
                                    num:marker_num, trimmed_label:the_m.marker_name.trimToPix(100),
                                    font_size:'font-size:0.95em;'
                        };
                        marker_num += 1;

                        var print_tpl = new Ext.Template(print_tpl_str);
                        var marker_html = print_tpl.applyTemplate(options);
                        html += marker_html;
                    }
                });
                margin = 21;
            }
            // If the img_url begins with ? it contains only parameters for ncfetch, so prepend ncfetch URL
            // This is a way to provide reliable URL resolution for embedding. SV-1760
            var img_url = from_cgi.img_url;
            if (img_url && img_url.charAt(0) == '?') {
                img_url = app.m_CGIs.NetCache + img_url;
            }
            html += '<img aria-label="icon" style="border: 1px solid;margin-top:'+margin+'px;" src="' + img_url + '"></div>';
            this.generator.document.write(html);


        },this);

        if (this.text_data) {
            this.generator.document.write('<br><br>Sequence View<br>');
            var top_nums = app.m_TextView.genTopNumbers();
            var html = SeqView.TextView.getViewTmpl().apply({border:1, width:img_width, left_nums:this.text_data.starts, top_nums:top_nums, id: app.m_Idx, sequence:this.text_data.sequence});
            this.generator.document.write(html);
            app.m_TextView.m_SequenceCols = this.col_tmp; // restore back
        }

        this.generator.document.write('</body></html>');
        this.generator.document.close();

    }

};
function globStringToRegexStr(str) {
    return preg_quote(str).replace(/\\\*/g, '.*').replace(/\\\?/g, '.');
}

function globStringToRegex(str) {
    return new RegExp(preg_quote(str).replace(/\\\*/g, '.*').replace(/\\\?/g, '.'), 'i');
}
function preg_quote (str, delimiter) {
    // http://kevin.vanzonneveld.net
    // +   original by: booeyOH
    // +   improved by: Ates Goral (http://magnetiq.com)
    // +   improved by: Kevin van Zonneveld (http://kevin.vanzonneveld.net)
    // +   bugfixed by: Onno Marsman
    // +   improved by: Brett Zamir (http://brett-zamir.me)
    // *     example 1: preg_quote("$40");
    // *     returns 1: '\$40'
    // *     example 2: preg_quote("*RRRING* Hello?");
    // *     returns 2: '\*RRRING\* Hello\?'
    // *     example 3: preg_quote("\\.+*?[^]$(){}=!<>|:");
    // *     returns 3: '\\\.\+\*\?\[\^\]\$\(\)\{\}\=\!\<\>\|\:'
    return (str + '').replace(new RegExp('[.\\\\+*?\\[\\^\\]$(){}=!<>|:\\' + (delimiter || '') + '-]', 'g'), '\\$&');
}

/** Polyfills for IE */
var $jscomp=$jscomp||{};
$jscomp.scope={};$jscomp.ASSUME_ES5=!1;$jscomp.ASSUME_NO_NATIVE_MAP=!1;$jscomp.ASSUME_NO_NATIVE_SET=!1;
$jscomp.defineProperty=$jscomp.ASSUME_ES5||"function"==typeof Object.defineProperties?Object.defineProperty:function(a,b,c){a!=Array.prototype&&a!=Object.prototype&&(a[b]=c.value)};
$jscomp.getGlobal=function(a){return"undefined"!=typeof window&&window===a?a:"undefined"!=typeof global&&null!=global?global:a};
$jscomp.global=$jscomp.getGlobal(this);$jscomp.SYMBOL_PREFIX="jscomp_symbol_";
$jscomp.initSymbol=function(){$jscomp.initSymbol=function(){};$jscomp.global.Symbol||($jscomp.global.Symbol=$jscomp.Symbol)};
$jscomp.Symbol=function(){var a=0;return function(b){return $jscomp.SYMBOL_PREFIX+(b||"")+a++}}();
$jscomp.initSymbolIterator=function(){$jscomp.initSymbol();var a=$jscomp.global.Symbol.iterator;a||(a=$jscomp.global.Symbol.iterator=$jscomp.global.Symbol("iterator"));"function"!=typeof Array.prototype[a]&&$jscomp.defineProperty(Array.prototype,a,{configurable:!0,writable:!0,value:function(){return $jscomp.arrayIterator(this)}});$jscomp.initSymbolIterator=function(){}};
$jscomp.arrayIterator=function(a){var b=0;return $jscomp.iteratorPrototype(function(){return b<a.length?{done:!1,value:a[b++]}:{done:!0}})};
$jscomp.iteratorPrototype=function(a){$jscomp.initSymbolIterator();a={next:a};a[$jscomp.global.Symbol.iterator]=function(){return this};return a};
$jscomp.iteratorFromArray=function(a,b){$jscomp.initSymbolIterator();a instanceof String&&(a+="");var c=0,d={next:function(){if(c<a.length){var e=c++;return{value:b(e,a[e]),done:!1}}d.next=function(){return{done:!0,value:void 0}};return d.next()}};d[Symbol.iterator]=function(){return d};return d};
$jscomp.polyfill=function(a,b,c,d){if(b){c=$jscomp.global;a=a.split(".");for(d=0;d<a.length-1;d++){var e=a[d];e in c||(c[e]={});c=c[e]}a=a[a.length-1];d=c[a];b=b(d);b!=d&&null!=b&&$jscomp.defineProperty(c,a,{configurable:!0,writable:!0,value:b})}};
$jscomp.polyfill("Array.prototype.values",function(a){return a?a:function(){return $jscomp.iteratorFromArray(this,function(a,c){return c})}},"es6","es3");
$jscomp.polyfill("Array.from",function(a){return a?a:function(a,c,d){$jscomp.initSymbolIterator();c=null!=c?c:function(a){return a};var b=[],f=a[Symbol.iterator];if("function"==typeof f)for(a=f.call(a);!(f=a.next()).done;)b.push(c.call(d,f.value));else{f=a.length;for(var g=0;g<f;g++)b.push(c.call(d,a[g]))}return b}},"es6","es3");
$jscomp.makeIterator=function(a){$jscomp.initSymbolIterator();var b=a[Symbol.iterator];return b?b.call(a):$jscomp.arrayIterator(a)};
$jscomp.FORCE_POLYFILL_PROMISE=!1;
$jscomp.polyfill("Promise",function(a){function b(){this.batch_=null}function c(a){return a instanceof e?a:new e(function(b,g){b(a)})}if(a&&!$jscomp.FORCE_POLYFILL_PROMISE)return a;b.prototype.asyncExecute=function(a){null==this.batch_&&(this.batch_=[],this.asyncExecuteBatch_());this.batch_.push(a);return this};b.prototype.asyncExecuteBatch_=function(){var a=this;this.asyncExecuteFunction(function(){a.executeBatch_()})};var d=$jscomp.global.setTimeout;b.prototype.asyncExecuteFunction=function(a){d(a,
0)};b.prototype.executeBatch_=function(){for(;this.batch_&&this.batch_.length;){var a=this.batch_;this.batch_=[];for(var b=0;b<a.length;++b){var c=a[b];delete a[b];try{c()}catch(l){this.asyncThrow_(l)}}}this.batch_=null};b.prototype.asyncThrow_=function(a){this.asyncExecuteFunction(function(){throw a;})};var e=function(a){this.state_=0;this.result_=void 0;this.onSettledCallbacks_=[];var b=this.createResolveAndReject_();try{a(b.resolve,b.reject)}catch(k){b.reject(k)}};e.prototype.createResolveAndReject_=
function(){function a(a){return function(d){c||(c=!0,a.call(b,d))}}var b=this,c=!1;return{resolve:a(this.resolveTo_),reject:a(this.reject_)}};e.prototype.resolveTo_=function(a){if(a===this)this.reject_(new TypeError("A Promise cannot resolve to itself"));else if(a instanceof e)this.settleSameAsPromise_(a);else{a:switch(typeof a){case "object":var b=null!=a;break a;case "function":b=!0;break a;default:b=!1}b?this.resolveToNonPromiseObj_(a):this.fulfill_(a)}};e.prototype.resolveToNonPromiseObj_=function(a){var b=
void 0;try{b=a.then}catch(k){this.reject_(k);return}"function"==typeof b?this.settleSameAsThenable_(b,a):this.fulfill_(a)};e.prototype.reject_=function(a){this.settle_(2,a)};e.prototype.fulfill_=function(a){this.settle_(1,a)};e.prototype.settle_=function(a,b){if(0!=this.state_)throw Error("Cannot settle("+a+", "+b|"): Promise already settled in state"+this.state_);this.state_=a;this.result_=b;this.executeOnSettledCallbacks_()};e.prototype.executeOnSettledCallbacks_=function(){if(null!=this.onSettledCallbacks_){for(var a=
this.onSettledCallbacks_,b=0;b<a.length;++b)a[b].call(),a[b]=null;this.onSettledCallbacks_=null}};var f=new b;e.prototype.settleSameAsPromise_=function(a){var b=this.createResolveAndReject_();a.callWhenSettled_(b.resolve,b.reject)};e.prototype.settleSameAsThenable_=function(a,b){var c=this.createResolveAndReject_();try{a.call(b,c.resolve,c.reject)}catch(l){c.reject(l)}};e.prototype.then=function(a,b){function c(a,b){return"function"==typeof a?function(b){try{d(a(b))}catch(m){f(m)}}:b}var d,f,h=new e(function(a,
b){d=a;f=b});this.callWhenSettled_(c(a,d),c(b,f));return h};e.prototype.catch=function(a){return this.then(void 0,a)};e.prototype.callWhenSettled_=function(a,b){function c(){switch(d.state_){case 1:a(d.result_);break;case 2:b(d.result_);break;default:throw Error("Unexpected state: "+d.state_);}}var d=this;null==this.onSettledCallbacks_?f.asyncExecute(c):this.onSettledCallbacks_.push(function(){f.asyncExecute(c)})};e.resolve=c;e.reject=function(a){return new e(function(b,c){c(a)})};e.race=function(a){return new e(function(b,
d){for(var e=$jscomp.makeIterator(a),f=e.next();!f.done;f=e.next())c(f.value).callWhenSettled_(b,d)})};e.all=function(a){var b=$jscomp.makeIterator(a),d=b.next();return d.done?c([]):new e(function(a,e){function f(b){return function(c){h[b]=c;g--;0==g&&a(h)}}var h=[],g=0;do h.push(void 0),g++,c(d.value).callWhenSettled_(f(h.length-1),e),d=b.next();while(!d.done)})};return e},"es6","es3");
$jscomp.owns=function(a,b){return Object.prototype.hasOwnProperty.call(a,b)};
$jscomp.polyfill("Object.assign",function(a){return a?a:function(a,c){for(var b=1;b<arguments.length;b++){var e=arguments[b];if(e)for(var f in e)$jscomp.owns(e,f)&&(a[f]=e[f])}return a}},"es6","es3");
$jscomp.polyfill("Array.prototype.keys",function(a){return a?a:function(){return $jscomp.iteratorFromArray(this,function(a){return a})}},"es6","es3");
$jscomp.polyfill("Array.from", function(a) {
    return a ? a : function(a, p, d) {
        $jscomp.initSymbolIterator();
        p = null != p ? p : function(a) {
            return a
        };
        var c = [],
            u = a[Symbol.iterator];
        if ("function" == typeof u)
            for (a = u.call(a); !(u = a.next()).done;)
                c.push(p.call(d, u.value));
        else {
            u = a.length;
            for (var z = 0; z < u; z++)
                c.push(p.call(d, a[z]))
        }
        return c
    }
}, "es6", "es3");

/** End of polyfill */

Ext.apply(SeqView.App, {

    /** API to download SV data
        @function downloadTracks
        @param  {Object}    request                     - download request
        @prop   {string}    id                          - sequence id
        @prop   {string}    format                      - data format (default: 'gff3')
        @prop   {Array}     tracks                      - tracks
        @prop   {Array}     ranges                      - ranges, array of [start, stop], 0-based
        @prop   {string?}   affinity_token              - affinity token, optional
        @prop   {Array?}    columns                     - columns   (TBD)
        @prop   {boolean}   transcr_prot_features       - include transcript & protein features
        @prop   {string?}   filename                    - file name, optional
        @prop   {Number?}   wait                        - debug: wait 'wait' ms before passing the response
        @param  {function?} callback                    - callback function, optional
        @type   {Promise}   loader                      - returns Promise with the following properties
        @prop   {string}    loader.url                  - current URL
        @prop   {boolean}   loader.cancelled            - whether request is cancelled
        @prop   {string}    loader.job                  - current job id
        @prop   {string}    loader.uuid
        @prop   {string}    loader.filename             - file name
        @prop   {function}  loader.cancel               - call to cancel request */
    downloadTracks: function(request, callback) {
        console.assert(request != null, "Parameter 'request' is required");
        var format = request.format || 'gff3',
            id = request.id,
            tracks = request.tracks || [],
            type = tracks.length > 0 ? tracks[0].key : null,
            ranges = request.ranges || [],
            filename = request.filename,
            assm_context = request.assm_context,
            wait = request.wait,        // for debug, wait ms
            target,
            cancelled,
            _cancel,
            job,
            _fileSize,      // file size
            url;

        // Get SV path
        var sv = Array.from(document.querySelectorAll("script")).filter(function(d){return d.src.match(/\/sviewer\.*\.js/)}).map(function(d){return d.src})[0],
            a = sv ? sv.match(/(https?:\/\/[a-z0-9.]+(?::[0-9]+)?\/.*)\/js\/sviewer\.js/i) : null,
            path = a ? a[1] : "/projects/sviewer";

        var t0 = Date.now(),
            getDelay = function(p){
                var dt = Date.now() - t0;
                p = p || 0;
                return Math.min(10000, Math.max(1000, 0.5 * (1 - p) * dt)); // it remains about (1-p)*dt*(1/p)
            },
            _wait = function(data, sec){
                return new Promise(function(resolve){
                    setTimeout(function(){
                        resolve(data);
                    }, sec);
                });
            },
            _callback = function(data){
                if(typeof callback == 'function'){
                    try{
                        callback(Object.assign({}, data, {target: target}));
                    }catch(e){
                        console.error(e);
                    }
                }
            };

        var _call = function(data){
                return new Promise(function(resolve, reject){
                var tm,
                    cancel = function(){
                        cancelled = true;
                        _cancel = null;
                        if(job){
                            if(tm){
                                clearTimeout(tm);
                                tm = null;
                            }
                            _fetch({
                                job_key: job,
                                'Cancel': 'Cancel'
                            });
                        }
                        reject(new Error("Cancelled"));
                    };

                var _fetch = function(data){
                    fetch(url, {
                        method: 'POST',
                        mode: 'cors',
                        cache: 'no-cache',                  // *default, no-cache, reload, force-cache, only-if-cached
                        credentials: 'same-origin',         // include, *same-origin, omit
                        headers: {
                          //'Content-Type': 'application/json'
                          'Content-Type': 'application/x-www-form-urlencoded'
                        },
                        redirect: 'follow',                 // manual, *follow, error
                        referrer: 'client',            // no-referrer, *client
                        body: Object.keys(data).filter(function(k){return data[k] != null}).map(function(k){return k + '=' + encodeURIComponent(data[k])}).join('&')   // URL encoded
                        // NCBI cannot handle JSON
                        //body: JSON.stringify(data)        // body data type must match "Content-Type" header
                    })
                    .then(function(response){
                        if(!response.ok){
                            throw new Error('HTTP Error  ' + response.status);
                        }
                        return response.json();
                    })
                    .then(function(data){
                        // debug
                        return wait > 0 ? _wait(data, wait) : data;
                    })
                    .then(function(data){
                        tm = null;
                        switch(data.job_status || 'done'){
                          case 'submitted':
                          case 'pending':
                          case 'running':
                            if(!cancelled){
                                _callback(data);

                                job = data.job_id;
                                tm = setTimeout(function(){
                                    tm = null;
                                    _fetch({
                                        job_key: data.job_id
                                    });
                                }, getDelay());
                            }
                            break;

                          case 'cancelled':
                          case 'canceled':
                            if(!_cancelled){    // job cancelled by other means
                                cancelled = true;
                                _cancel = null;
                                reject(new Error("Cancelled"));
                            }
                            break;

                          case 'done':
                            _cancel = null;     // reset cancel
                            resolve(data);
                            break;

                          case 'failed':
                            throw new Error(data.error_message || "Unknown error");
                            break;

                          default:
                            throw new Error('Invalid job status' +  data.job_status);
                            break;
                        }
                    })
                    .catch(function(e){
                        _cancel = null;
                        reject(e);
                    });
                };

                _cancel = cancel;

                _fetch(data);
            });
        };

        var create_UUID = function(){
            var t = t0;
            return "xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx".replace(/[xy]/g, function(c){
                var r = (t + Math.random()*16)%16 | 0;
                t = Math.floor(t/16);
                return (c=='x' ? r :(r&0x3|0x8)).toString(16);
            });
            return uuid;
        };

        var encodeTracks = function(tracks){return SeqView.TM.tracksArrayToString(tracks || [])};

        var _bool = function(v){return typeof v == 'string' ? v == 'true' || v == '1' : Boolean(v)};

        var _create = function(){

            if (tracks[0].data_provider.indexOf('_rmt' > 0))
                url = path + '/'+tracks[0].data_provider+'.cgi';
            else
                url = path + '/seqgraphic.cgi';

            target = 'asn1';
            job = null;
            return _call({
                  target: "asn1",
                  id: id,
                  ranges: ranges.map(function(r){return r[0] + '-' + r[1]}).join(','),
                  affinity_token: request.affinity_token,
                  tracks: encodeTracks(tracks),
                  assm_context :assm_context
            })
            .then(function(data){
                if(_bool(data.success)){
                    url = path + '/sequence.cgi';
                    target = 'sequence';
                    job = null;
                    return _call({
                        blob_type : data.blob_type ,
                        blob_key : data.blob_key ,
                        ranges: ranges.map(function(r){return r[0] + '-' + r[1]}).join(','),
                        columns: request.columns,
                        transcr_prot_features: request.transcr_prot_features,
                        track_type: type,
                        id:id,
                        format: format
                    });
                }else{
                    return Promise.reject(new Error(data.msg));
                }
            })
            .then(function(data){
                _fileSize = data._file_size;
                return _bool(data.success)
                    ? Promise.resolve(data)
                    : Promise.reject(new Error(data.msg));
            });
        };

        var _download = new Promise(function(resolve, reject){
            _create()
            .then(function(data){
                var track = tracks[0],
                    title = track ? track.display_name : null,
                    n = 15,
                    f = false,
                    //f = navigator.platform != 'MacIntel', // do not truncate on Mac
                    _trunc = function(s){return f && s && s.length > 2*n + 1 ? s.substr(0, n) + '_' + s.substr(s.length - n) : s};
                var pref_ext = data.file_ext && data.file_ext.length ? format + '.' + data.file_ext : format;
                filename = filename || (title ? _trunc(title) + '.' + pref_ext : 'table.' + pref_ext);
                filename = filename.replace(/\,/g, ' ');
                var url = path + '/ncfetch.cgi' + data.file_url;

                setTimeout(function(){
                    var a = document.createElement("a");
                    a.href = url + '&filename=' + filename;
                    a.download = filename; // ignored with cross-origin

                    document.body.appendChild(a);
                    a.click();

                    document.body.removeChild(a);
                }, 100);
                resolve(data);
            })
            .catch(function(e){
                reject(e);
            });
        });

        Object.defineProperties(_download, {
            cancel: {
                value: function(){
                    cancelled = true;
                    if(_cancel){
                        _cancel();
                    }
                }
            },
            cancelled: {
                get: function(){
                    return Boolean(cancelled);
                }
            },
            url: {
                get: function(){
                    return url;
                }
            },
            fileSize: {
                get: function(){
                    return _fileSize;
                }
            },
            filename: {
                get: function(){
                    return filename;
                }
            },
            job: {
                get: function(){
                    return job;
                }
            },
            uuid: {
                value: create_UUID()
            }
        });

        return _download;
    }

});
}
