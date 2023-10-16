
// This script collapses and expands a section by clicking on it.
// This is used in the "optional parameters" section of the umIT's functions documentation.
var coll = document.getElementsByClassName("mycollapsible");
var i;

for (i = 0; i < coll.length; i++) {
  coll[i].addEventListener("click", function() {
    this.classList.toggle("active");
    var content = this.nextElementSibling;
    var icon = document.getElement
    if (content.style.display === "block") {
      content.style.display = "none";
      console.log("closed")
      this.getElementsByClassName("iconminus")[0].className = "iconplus";
    } else {
      content.style.display = "block";
      console.log("open")
      this.getElementsByClassName("iconplus")[0].className = "iconminus";
    }    
  });
}
