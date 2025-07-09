{
  description = "Flake to run new_scf";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = {
    self,
    nixpkgs,
    flake-utils,
    ...
  }:
    flake-utils.lib.eachDefaultSystem (
      system: let
        pkgs = import nixpkgs {
          inherit system;
        };
      in {
        packages.default = pkgs.stdenv.mkDerivation {
          name = "new_scf";

          src = builtins.fetchGit {
            url = "https://github.com/prajvalk/new_scf.git";
            rev = "c6b24f512a440a66fb1f146c5d0e1d5531aba38c";
            submodules = true;
          };

          nativeBuildInputs = with pkgs; [
            gcc
            gnumake
            cmake
            openblas
          ];

          configurePhase = ''
            mkdir build
            cd build
            cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++
            cd ..
          '';

          buildPhase = ''
            cd build
            make -j8
          '';
        };

        apps.default = flake-utils.lib.mkApp {
          drv = self.packages.${system}.default;
        };
      }
    );
}
