const MAX = 0;
const MEDIAN = 1;
const COMPARISON = 2;

function disableElement(id) {
    document.getElementById(id).disabled = true;
}

function enableElement(id) {
    document.getElementById(id).disabled = false;
}

function checkboxClicked(checkbox_name) {
    let checkbox = document.getElementsByName(checkbox_name)[0];
    if (checkbox.checked)
        enableElement("t0");
    else
        disableElement("t0");
}

function changeMode(mode_name) {
    const mode = document.getElementsByName(mode_name)[0].selectedIndex;
    let input = document.getElementById("t0");

    switch (mode) {
        case MAX:
            disableElement("g2");
            //disableElement("t0");
            input.setAttribute("value", "500");
            break;
        case MEDIAN:
            disableElement("g2");
            //enableElement("t0");
            input.setAttribute("value", "2.0");
            break;
        case COMPARISON:
            enableElement("g2");
            //enableElement("t0");
            input.setAttribute("value", "2.0");
            break;
    }
}

function checkAGIs() {
    // Note regId is global. It should work.
    const primary = document.getElementById("g1").value;
    const secondary = document.getElementById("g2").value;
    const mode = document.getElementsByName("mode")[0].selectedIndex;
    if (primary === "") {
        alert("The Primary AGI Field Cannot Be Empty!");
        return false;
    } else if (primary.match(regId) === null) {
        alert("Invalid Primary AGI: Not in Correct Format!");
        return false;
    } else if (mode === COMPARISON) {
        if (secondary === "") {
            alert("The Secondary AGI Field Cannot Be Empty in Compare Mode!");
            return false;
        } else if (secondary.match(regId) == null) {
            alert("Invalid Secondary AGI: Not in Correct Format!");
            return false;
        } else if (secondary === primary) {
            alert("Primary and Secondary genes should be different!");
            return false;
        }
    }
    return true;
}

function openChart(element) {
    const layer = document.getElementById(element);
    layer.style.visibility = 'visible';
}

function closeChart(element) {
    const layer = document.getElementById(element);
    layer.style.visibility = 'hidden';
}

function resizeIframe(iframeId, iframe) {
    const target = document.getElementById(iframeId);
    target.height = iframe.document.body.offsetHeight + 25; // Firefox
    target.style.height = iframe.document.body.scrollHeight + 25 // Opera and IE
}
