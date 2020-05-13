function unitConverter() {

  var hbar = 1.054571817e-34;
  var c = 299792458.0;
  var eV = 1.602176634e-19;

  var length = hbar*c/eV;
  var time = hbar/eV;

  var inputValue = document.getElementById("inputValue").value;
  var unit = document.getElementById("inputUnit").value;

  if (unit == "Hz") {
    document.getElementById("outputValue").innerHTML = inputValue * time;
  } else if (unit == "m") {
    document.getElementById("outputValue").innerHTML = inputValue / length;
  } else if (unit == "eV") {
    document.getElementById("outputValue").innerHTML = inputValue;
  }
}
function unitConverterFreq() {

  var hbar = 1.054571817e-34;
  var c = 299792458.0;
  var eV = 1.602176634e-19;

  var length = hbar*c/eV;
  var time = hbar/eV;

  var inputValue = document.getElementById("inputValueFreq").value;
  var unit = document.getElementById("inputUnitFreq").value;

  if (unit == "Hz") {
    document.getElementById("outputValueFreq").innerHTML = inputValue * time;
  } else if (unit == "m") {
    document.getElementById("outputValueFreq").innerHTML = inputValue / length;
  } else if (unit == "eV") {
    document.getElementById("outputValueFreq").innerHTML = inputValue;
  }
}
function unitConverterPol() {

  var hbar = 1.054571817e-34;
  var c = 299792458.0;
  var eV = 1.602176634e-19;

  var length = hbar*c/eV;
  var latom = 5.292e-11;
  var lnat  = 1.97e-7;
  var  langstrom = 1e-10;

  var inputValue = document.getElementById("inputValuePol").value;
  var unit = document.getElementById("inputUnitPol").value;

  if (unit == "ang3") {
    document.getElementById("outputValuePol").innerHTML = inputValue * (langstrom/lnat)**3;
  } else if (unit == "au") {
    document.getElementById("outputValuePol").innerHTML = inputValue * (latom/lnat)**3;
  } else if (unit == "eV-3") {
    document.getElementById("outputValuePol").innerHTML = inputValue;
  }
}
function unitConverterLen() {

  var hbar = 1.054571817e-34;
  var c = 299792458.0;
  var eV = 1.602176634e-19;

  var length = hbar*c/eV;
  var time = hbar/eV;

  var inputValue = document.getElementById("inputValueLen").value;
  var unit = document.getElementById("inputUnitLen").value;

  if (unit == "m") {
    document.getElementById("outputValueLen").innerHTML = inputValue / length;
  } else if (unit == "eV-1") {
    document.getElementById("outputValueLen").innerHTML = inputValue;
  }
}
function unitConverterVel() {

  var c = 299792458.0;

  var inputValue = document.getElementById("inputValueVel").value;
  var unit = document.getElementById("inputUnitVel").value;

  if (unit == "m/s") {
    document.getElementById("outputValueVel").innerHTML = inputValue / c;
  } else if (unit == "c") {
    document.getElementById("outputValueVel").innerHTML = inputValue;
  }
}
function unitConverterTemp() {

  var nattemp = 1.16e4;

  var inputValue = document.getElementById("inputValueTemp").value;
  var unit = document.getElementById("inputUnitTemp").value;

  if (unit == "1/K") {
    document.getElementById("outputValueTemp").innerHTML = inputValue * nattemp;
  } else if (unit == "eV-1") {
    document.getElementById("outputValueTemp").innerHTML = inputValue;
  }
}
