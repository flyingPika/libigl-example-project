/*
#include <igl/readOFF.h>
//#undef IGL_STATIC_LIBRARY
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <igl/opengl/glfw/Viewer.h>

#include <Eigen/Core>
#include <iostream>


Eigen::MatrixXd VA, VB, VC;
Eigen::VectorXi J, I;
Eigen::MatrixXi FA, FB, FC;
igl::MeshBooleanType boolean_type(
    igl::MESH_BOOLEAN_TYPE_UNION);

const char* MESH_BOOLEAN_TYPE_NAMES[] =
{
  "Union",
  "Intersect",
  "Minus",
  "XOR",
  "Resolve",
};

void update(igl::opengl::glfw::Viewer& viewer)
{
    igl::copyleft::cgal::mesh_boolean(VA, FA, VB, FB, boolean_type, VC, FC, J);
    Eigen::MatrixXd C(FC.rows(), 3);
    for (size_t f = 0; f < C.rows(); f++)
    {
        if (J(f) < FA.rows())
        {
            C.row(f) = Eigen::RowVector3d(1, 0, 0);
        }
        else
        {
            C.row(f) = Eigen::RowVector3d(0, 1, 0);
        }
    }
    viewer.data().clear();
    viewer.data().set_mesh(VC, FC);
    viewer.data().set_colors(C);
    std::cout << "A " << MESH_BOOLEAN_TYPE_NAMES[boolean_type] << " B." << std::endl;
}

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int mods)
{
    switch (key)
    {
    default:
        return false;
    case '.':
        boolean_type =
            static_cast<igl::MeshBooleanType>(
                (boolean_type + 1) % igl::NUM_MESH_BOOLEAN_TYPES);
        break;
    case ',':
        boolean_type =
            static_cast<igl::MeshBooleanType>(
                (boolean_type + igl::NUM_MESH_BOOLEAN_TYPES - 1) %
                igl::NUM_MESH_BOOLEAN_TYPES);
        break;
    case '[':
        viewer.core().camera_dnear -= 0.1;
        return true;
    case ']':
        viewer.core().camera_dnear += 0.1;
        return true;
    }
    update(viewer);
    return true;
}

int main(int argc, char* argv[])
{
    using namespace Eigen;
    using namespace std;
    igl::readOFF("D:\\myDocuments\\projects\\github\\libigl-example-project\\Models\\fertility.off", VA, FA);
    igl::readOFF("D:\\myDocuments\\projects\\github\\libigl-example-project\\Models\\kitten.off", VB, FB);
    // Plot the mesh with pseudocolors
    igl::opengl::glfw::Viewer viewer;

    // Initialize
    update(viewer);

    viewer.data().show_lines = true;
    viewer.callback_key_down = &key_down;
    viewer.core().camera_dnear = 3.9;
    cout <<
        "Press '.' to switch to next boolean operation type." << endl <<
        "Press ',' to switch to previous boolean operation type." << endl <<
        "Press ']' to push near cutting plane away from camera." << endl <<
        "Press '[' to pull near cutting plane closer to camera." << endl <<
        "Hint: investigate _inside_ the model to see orientation changes." << endl;
    viewer.launch();
}
*/

#include <igl/avg_edge_length.h>
#include <igl/cotmatrix.h>
#include <igl/invert_diag.h>
#include <igl/massmatrix.h>
#include <igl/parula.h>
#include <igl/per_corner_normals.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/principal_curvature.h>
#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

int main(int argc, char* argv[])
{
    using namespace Eigen;
    std::string filename = "C:/Users/zhaoj/Desktop/vase100K.off";
    if (argc > 1)
    {
        filename = argv[1];
    }
    // Load a mesh in OFF format
    igl::read_triangle_mesh(filename, V, F);

    // Alternative discrete mean curvature
    MatrixXd HN;
    SparseMatrix<double> L, M, Minv;
    igl::cotmatrix(V, F, L);
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
    igl::invert_diag(M, Minv);
    // Laplace-Beltrami of position
    HN = -Minv * (L * V);
    // Extract magnitude as mean curvature
    VectorXd H = HN.rowwise().norm();

    // Compute curvature directions via quadric fitting
    MatrixXd PD1, PD2;
    VectorXd PV1, PV2;
    igl::principal_curvature(V, F, PD1, PD2, PV1, PV2);
    // mean curvature
    H = 0.5 * (PV1 + PV2);

    igl::opengl::glfw::Viewer viewer;
    viewer.data().set_mesh(V, F);

    viewer.data().set_data(H);

    // Average edge length for sizing
    const double avg = igl::avg_edge_length(V, F);

    // Draw a red segment parallel to the maximal curvature direction
    // const RowVector3d red(0.8, 0.2, 0.2), blue(0.2, 0.2, 0.8);
    // viewer.data().add_edges(V + PD1 * avg, V - PD1 * avg, red);

    // Draw a blue segment parallel to the minimal curvature direction
    // viewer.data().add_edges(V + PD2 * avg, V - PD2 * avg, blue);

    // Hide wireframe
    viewer.data().show_lines = false;

    viewer.launch();
}