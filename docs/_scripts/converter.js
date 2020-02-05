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
  }
}
