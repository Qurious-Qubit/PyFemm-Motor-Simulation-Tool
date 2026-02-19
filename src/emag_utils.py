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



'''
    def create_operating_point(self, op_name, current, angle, speed=None):
        """Create a subdirectory for a specific operating point"""
        op_dir = self.operating_points_dir / op_name
        op_dir.mkdir(exist_ok=True)
        
        # Save operating point parameters
        op_config = {
            'name': op_name,
            'current': current,
            'angle': angle,
            'speed': speed if speed else (120 * self.frequency) / (2 * self.polepairs),
            'timestamp': datetime.now().isoformat()
        }
        
        # You can save this as JSON or any other format
        config_file = op_dir / "config.json"
        # json.dump(op_config, open(config_file, 'w'), indent=2)
        
        print(f"✓ Operating point created: {op_dir}")
        return op_dir
    
    def get_simulation_paths(self):
        """Return all important paths for easy access"""
        return {
            'base_path': self.base_path,
            'femfiles_dir': self.femfiles_dir,
            'simulations_dir': self.simulations_dir,
            'simulation_dir': self.simulation_dir,
            'femfile_copy_dir': self.femfile_copy_dir,
            'operating_points_dir': self.operating_points_dir,
            'fem_file_path': self.fem_file_path
        }
    
    def display_structure(self):
        """Display the directory structure"""
        print("\n" + "="*50)
        print("MOTOR SIMULATION DIRECTORY STRUCTURE")
        print("="*50)
        paths = self.get_simulation_paths()
        for name, path in paths.items():
            exists = "✓ EXISTS" if path.exists() else "⨯ MISSING"
            print(f"{name:25} : {path} {exists}")


def SimDC():

    return 1

'''
       

path = r"C:\D_Drive\Technical\GitHub\PyFemm-Motor-Simulation-Tool\test"
femmfilename = "example"
polepairs = 4
a = Motor(path = path, femmfilename = femmfilename, polepairs = polepairs)

