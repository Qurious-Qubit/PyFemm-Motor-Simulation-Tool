#This library contains utility functions for motor simulations using FEMM with post-processing in Python.
#This is specific for 3Phase motors and generalisation is future work.

import femm
import os
import shutil
from pathlib import Path

#assumed that the required inputs are in a class object format


class Motor:
    def __init__(self, **kwargs):
        # REQUIRED parameters
        required_params = ['path', 'femmfilename']
        
        # Check required parameters first
        for param in required_params:
            if param not in kwargs:
                raise ValueError(f"Missing required parameter: {param}")
        
        # Set ALL kwargs as direct attributes
        for key, value in kwargs.items():
            setattr(self, key, value)
        
        # Convert path to Path object
        self.path = Path(self.path)
        
        # Initialize directory structure
        self._setup_directory_structure()
    
    def _setup_directory_structure(self):
        """Create the required directory structure"""
        # Define all paths
        self.femfiles_dir = self.path / "femfiles"
        self.simulations_dir = self.path / "simulations"
        self.simulation_dir = self.simulations_dir / self.femmfilename
        self.femfile_copy_dir = self.simulation_dir / "femfile"
        self.operating_points_dir = self.simulation_dir / "OperatingPoints"
        
        # Create directories if they don't exist
        directories = [
            self.femfiles_dir,
            self.simulations_dir, 
            self.simulation_dir,
            self.femfile_copy_dir,
            self.operating_points_dir
        ]
        
        for directory in directories:
            directory.mkdir(parents=True, exist_ok=True)
            print(f"✓ Directory ensured: {directory}")
 
        # Copy FEM file if it exists in femfiles directory
        self._copy_fem_file()
        
    def _copy_fem_file(self):
        """Copy the FEM file from femfiles to simulation directory"""
        source_fem_file = self.femfiles_dir / f"{self.femmfilename}.fem"
        destination_fem_file = self.femfile_copy_dir / f"{self.femmfilename}.fem"
        
        if source_fem_file.exists():
            shutil.copy2(source_fem_file, destination_fem_file)
            print(f"✓ FEM file copied: {source_fem_file} → {destination_fem_file}")
            self.fem_file_path = destination_fem_file
        else:
            raise FileNotFoundError(
                f"FEM file not found: {source_fem_file}\n"
                f"Please ensure the file exists in the 'femfiles' directory.\n"
                f"Expected path: {source_fem_file}"
            )

path = r"C:\D_Drive\Technical\GitHub\PyFemm-Motor-Simulation-Tool\test"
femmfilename = "example"
polepairs = 4
a = Motor(path = path, femmfilename = femmfilename, polepairs = polepairs)

