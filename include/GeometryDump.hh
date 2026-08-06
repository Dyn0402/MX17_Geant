#pragma once
// GeometryDump.hh
//
// Writes the constructed Geant4 geometry to a JSON file so the plotting
// scripts (scripts/model/plot_mx17_model.py) can render the TRUE geometry
// instead of a hand-maintained python mirror. Triggered by
// `mm_sim --dump-geometry <file>` after initialization; no beam is run.
//
// Format: a JSON array, one entry per placed physical volume (the world is
// skipped), positions in mm in world coordinates:
//   { "pv": ..., "lv": ..., "mat": ...,
//     "pos": [x, y, z],            // centre of the solid's local frame
//     "rot": [9 numbers],          // row-major, only when not identity
//     "solid": { "type": "box",  "hx":..,"hy":..,"hz":.. }
//           or { "type": "tubs", "rin":..,"rout":..,"hz":.. }
//           or { "type": "sub", "a": {...}, "b": {...}, "bpos": [x,y,z] }
//           or { "type": "other", "name": ..., "bb": [6 numbers] } }

#include <string>

class G4VPhysicalVolume;

void DumpGeometry(const G4VPhysicalVolume* world, const std::string& file);
