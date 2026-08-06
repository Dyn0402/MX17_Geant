// GeometryDump.cc — see GeometryDump.hh

#include "GeometryDump.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4DisplacedSolid.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4VisExtent.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>
#include <iomanip>
#include <sstream>

namespace {

void SolidJson(const G4VSolid* solid, std::ostream& os) {
    if (auto* box = dynamic_cast<const G4Box*>(solid)) {
        os << "{\"type\":\"box\",\"hx\":" << box->GetXHalfLength()/mm
           << ",\"hy\":" << box->GetYHalfLength()/mm
           << ",\"hz\":" << box->GetZHalfLength()/mm << "}";
        return;
    }
    if (auto* tub = dynamic_cast<const G4Tubs*>(solid)) {
        os << "{\"type\":\"tubs\",\"rin\":" << tub->GetInnerRadius()/mm
           << ",\"rout\":" << tub->GetOuterRadius()/mm
           << ",\"hz\":" << tub->GetZHalfLength()/mm << "}";
        return;
    }
    if (auto* sub = dynamic_cast<const G4SubtractionSolid*>(solid)) {
        const G4VSolid* a = sub->GetConstituentSolid(0);
        const G4VSolid* b = sub->GetConstituentSolid(1);
        G4ThreeVector bpos(0, 0, 0);
        if (auto* disp = dynamic_cast<const G4DisplacedSolid*>(b)) {
            bpos = disp->GetObjectTranslation();
            b = disp->GetConstituentMovedSolid();
        }
        os << "{\"type\":\"sub\",\"a\":";
        SolidJson(a, os);
        os << ",\"b\":";
        SolidJson(b, os);
        os << ",\"bpos\":[" << bpos.x()/mm << "," << bpos.y()/mm << ","
           << bpos.z()/mm << "]}";
        return;
    }
    const G4VisExtent e = solid->GetExtent();
    os << "{\"type\":\"other\",\"name\":\"" << solid->GetEntityType()
       << "\",\"bb\":[" << e.GetXmin()/mm << "," << e.GetXmax()/mm << ","
       << e.GetYmin()/mm << "," << e.GetYmax()/mm << ","
       << e.GetZmin()/mm << "," << e.GetZmax()/mm << "]}";
}

void Walk(const G4VPhysicalVolume* pv, const G4RotationMatrix& parentRot,
          const G4ThreeVector& parentPos, bool isWorld, std::ostream& os,
          bool& first) {
    // Placement transform: objectRotation is the ACTIVE rotation of the solid.
    const G4RotationMatrix rot = parentRot * pv->GetObjectRotationValue();
    const G4ThreeVector pos = parentPos + parentRot * pv->GetObjectTranslation();

    const G4LogicalVolume* lv = pv->GetLogicalVolume();
    if (!isWorld) {
        if (!first) os << ",\n";
        first = false;
        os << "{\"pv\":\"" << pv->GetName() << "\",\"lv\":\"" << lv->GetName()
           << "\",\"mat\":\"" << lv->GetMaterial()->GetName() << "\","
           << "\"pos\":[" << pos.x()/mm << "," << pos.y()/mm << ","
           << pos.z()/mm << "]";
        if (!rot.isIdentity()) {
            os << ",\"rot\":[" << rot.xx() << "," << rot.xy() << "," << rot.xz()
               << "," << rot.yx() << "," << rot.yy() << "," << rot.yz()
               << "," << rot.zx() << "," << rot.zy() << "," << rot.zz() << "]";
        }
        os << ",\"solid\":";
        SolidJson(lv->GetSolid(), os);
        os << "}";
    }
    for (size_t i = 0; i < lv->GetNoDaughters(); ++i)
        Walk(lv->GetDaughter(i), rot, pos, false, os, first);
}

}  // namespace

void DumpGeometry(const G4VPhysicalVolume* world, const std::string& file) {
    std::ofstream os(file);
    os << std::setprecision(10);
    os << "[\n";
    bool first = true;
    Walk(world, G4RotationMatrix(), G4ThreeVector(), true, os, first);
    os << "\n]\n";
    os.close();
    G4cout << "GeometryDump: wrote " << file << G4endl;
}
