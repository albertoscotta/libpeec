/**
 * \file
 * \author Scotta Alberto
 *
\verbatim
 |¯¯| |    |¯¯\ \  /
 |--| |    |--<  \/
 |  | |___ |__/  /   
\endverbatim
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "mesh.h"
#include "utils/lin_algebra.h"

void sort_face_nodes(int face_num, FACE* face);

int build_iedge_connectivity(MESH* mesh);

void get_face_geom(MESH* mesh);

void get_iedge_geom(MESH* mesh);

MESH* build_mesh(int node_num, int face_num, const double(* node)[3],
  const int(* face2node)[3]) {
  
  if(!(node && face2node)) {
    #ifdef DEBUG
    fprintf(stderr, "Error: build_mesh(), NULL pointer ");
    fprintf(stderr, "detected in arguments.\n");
    #endif
    return NULL;
  }

  MESH* mesh;
  mesh = NULL;

  do {
    if(!(mesh = malloc(sizeof(MESH)))) {
      #ifdef DEBUG
      fprintf(stderr, "Error: build_mesh(), malloc() failed.\n");
      #endif
      break;
    }
    mesh->node = NULL;
    mesh->iedge = NULL;
    mesh->face = NULL;

    mesh->node_num = node_num;
    mesh->face_num = face_num;

    mesh->node = malloc(mesh->node_num*sizeof(NODE));
    mesh->face = malloc(mesh->face_num*sizeof(FACE));
    if(!(mesh->node && mesh->face)){
      #ifdef DEBUG
      fprintf(stderr, "Error: build_mesh(), malloc() failed.\n");
      #endif  
      break;
    }

    memcpy(mesh->node, node, node_num*3*sizeof(double));
    #pragma omp parallel for
    for(int i = 0; i < face_num; i ++)
      memcpy(mesh->face[i].node, face2node[i], 3*sizeof(int));

    sort_face_nodes(mesh->face_num, mesh->face);

    if(build_iedge_connectivity(mesh))
      break;

    get_iedge_geom(mesh);
    get_face_geom(mesh);

    return mesh;
  } while(0);

  destroy_mesh(mesh);
  return NULL;
}

void destroy_mesh(MESH* mesh) {

  if(!mesh)
    return;

  if(mesh->node)
    free(mesh->node);
  if(mesh->iedge)
    free(mesh->iedge);
  if(mesh->face)
    free(mesh->face);
  free(mesh);

  return;
}

void sort_face_nodes(int face_num, FACE* face) {
  int i;
  int temp;

  #pragma omp parallel for private(temp)
  for(i = 0; i < face_num; i ++) {

    if(face[i].node[0] < face[i].node[1]) {
      if(face[i].node[2] < face[i].node[0]) { /* -> 201 */
        temp = face[i].node[0];
        face[i].node[0] = face[i].node[2];
        face[i].node[2] = face[i].node[1];
        face[i].node[1] = temp;
      } else {
        if(face[i].node[2] < face[i].node[1]) { /* -> 021 */
          temp = face[i].node[1];
          face[i].node[1] = face[i].node[2];
          face[i].node[2] = temp;
        }
      }
    } else {
      if(face[i].node[2] < face[i].node[1]) { /* -> 210 */
        temp = face[i].node[0];
        face[i].node[0] = face[i].node[2];
        face[i].node[2] = temp;
      } else {
        if(face[i].node[2] < face[i].node[0]) { /* -> 120 */
          temp = face[i].node[0];
          face[i].node[0] = face[i].node[1];
          face[i].node[1] = face[i].node[2];
          face[i].node[2] = temp;
        } else {              /* -> 102 */
          temp = face[i].node[0];
          face[i].node[0] = face[i].node[1];
          face[i].node[1] = temp;
        }
      }
    }
  }

  return;
}

int build_iedge_connectivity(MESH* mesh) {
  int bedge_num;
  int max_edge_num;
  int edge_nodes[2];
  int ex_bedge_local_index;
  int k;
  bool found;
  EDGE* edge;

  mesh->iedge_num = 0;
  bedge_num = 0;
  max_edge_num = 3*mesh->face_num;

  if(!(edge = malloc(max_edge_num*sizeof(EDGE)))) {
    #ifdef DEBUG  
    fprintf(stderr, "Error: build_iedge_connectivity(), malloc() failed,\n");
    #endif
    return 1;
  }

  for(int i = 0; i < mesh->face_num; i ++)
    for(int j = 0; j < 3; j ++) {
      mesh->face[i].iedge[j] = -1;

      edge_nodes[0] = mesh->face[i].node[(2-j)/2];
      edge_nodes[1] = mesh->face[i].node[(5-j)/2];

      {
        found = false;
        for(k = 0; k < bedge_num; k ++)
          if((edge_nodes[0] == edge[max_edge_num-1-k].node[0]) &&
             (edge_nodes[1] == edge[max_edge_num-1-k].node[1])) {
            found = true;
            break;
          }

        if(found) {
          edge[mesh->iedge_num].node[0] = edge_nodes[0];
          edge[mesh->iedge_num].node[1] = edge_nodes[1];

          edge[mesh->iedge_num].face[0] = edge[max_edge_num-1-k].face[0];
          edge[mesh->iedge_num].face[1] = i;

          mesh->face[i].iedge[j] = mesh->iedge_num;
          ex_bedge_local_index = edge[max_edge_num-1-k].face[1];
          mesh->face[edge[max_edge_num-1-k].face[0]].iedge[ex_bedge_local_index]
             = mesh->iedge_num;

          mesh->iedge_num ++;

          edge[max_edge_num-1-k].node[0] =
          edge[max_edge_num-bedge_num].node[0];
          edge[max_edge_num-1-k].node[1] =
          edge[max_edge_num-bedge_num].node[1];

          edge[max_edge_num-1-k].face[0] =
            edge[max_edge_num-bedge_num].face[0];
          edge[max_edge_num-1-k].face[1] =
            edge[max_edge_num-bedge_num].face[1];      

          bedge_num --;
        } else {
          bedge_num ++;

          edge[max_edge_num-bedge_num].node[0] = edge_nodes[0];
          edge[max_edge_num-bedge_num].node[1] = edge_nodes[1];

          edge[max_edge_num-bedge_num].face[0] = i;
          edge[max_edge_num-bedge_num].face[1] = j;
        }
      }
    }

  if(!(mesh->iedge = realloc(edge, mesh->iedge_num*sizeof(EDGE)))) {
    #ifdef DEBUG
    fprintf(stderr, "Error: build_iedge_connectivity(), realloc() failed.\n");
    #endif
    free(edge);
    return 1;
  }

  return 0;
}

void get_face_geom(MESH* mesh) {
  double edge_one_vector[3], edge_two_vector[3];
  double face_normal_vector[3];

  #pragma omp parallel for private(edge_one_vector, edge_two_vector, \
    face_normal_vector)
  for(int i = 0; i < mesh->face_num; i ++) {
    /* barycenter */
    vector_init(mesh->face[i].center);
    for(int j = 0; j < 3; j ++) /* loop over nodes */
      vector_add(mesh->face[i].center, mesh->face[i].center,
        (double *)&(mesh->node[mesh->face[i].node[j]]));
    vector_scale(mesh->face[i].center, 1.0/3.0);

    /* area */
    vector_sub(edge_one_vector,
      (double *)&(mesh->node[mesh->face[i].node[1]]),
      (double *)&(mesh->node[mesh->face[i].node[0]]));
    vector_sub(edge_two_vector,
      (double *)&(mesh->node[mesh->face[i].node[2]]),
      (double *)&(mesh->node[mesh->face[i].node[0]]));
    cross_product(face_normal_vector, edge_one_vector, edge_two_vector);
    mesh->face[i].area = sqrt(dot_product(face_normal_vector,
      face_normal_vector))/2.0;
  }

  return;
}

void get_iedge_geom(MESH* mesh) {
  double iedge_vector[3];

  #pragma omp parallel for private(iedge_vector)
  for(int i = 0; i < mesh->iedge_num; i ++) {
    /* midpoint */
    vector_add(mesh->iedge[i].midpoint,
      (double *)&(mesh->node[mesh->iedge[i].node[1]]),
      (double *)&(mesh->node[mesh->iedge[i].node[0]]));
    vector_scale(mesh->iedge[i].midpoint, 1.0/2.0);

    /* length */
    vector_sub(iedge_vector, (double *)&(mesh->node[mesh->iedge[i].node[1]]),
      (double *)&(mesh->node[mesh->iedge[i].node[0]]));
    mesh->iedge[i].length = sqrt(dot_product(iedge_vector, iedge_vector));
      
  }
  return;
}
