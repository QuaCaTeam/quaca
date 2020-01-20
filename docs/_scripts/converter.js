function frequencyConverter() {
  var valNum = document.getElementById("inputFrequency").value;
  var option = document.getElementById("inputFrequencyUnit").value;

  if (option == "eV") {
    document.getElementById("outputFrequency").innerHTML = valNum / 2;
  } else if (option == "nm") {
    document.getElementById("outputFrequency").innerHTML = valNum / 3;
  }
}
