{
  rPackages,
}:

rPackages.buildRPackage {
  pname = "mcclust.ext";
  version = "0.0.1";

  src = ./.;

  buildInputs = [
    rPackages.mcclust
  ];

  meta = {
    homepage = "https://github.com/sarawade/mcclust.ext";
  };
}
