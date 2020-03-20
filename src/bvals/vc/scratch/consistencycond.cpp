  if (pmb->block_size.nx3 == -1) {  // ignore
    _FinalizeVert3a();
    const int ng = NGHOST;

    // idx components for addressing block vertices and midpoints
    const int kcs = pmb->kvs, kcm = pmb->kvs+pmb->block_size.nx3/2, kce = pmb->kve;
    const int jcs = pmb->jvs, jcm = pmb->jvs+pmb->block_size.nx2/2, jce = pmb->jve;
    const int ics = pmb->ivs, icm = pmb->ivs+pmb->block_size.nx1/2, ice = pmb->ive;

    const int k_c[3] = {pmb->kvs, pmb->kvs+pmb->block_size.nx3/2, pmb->kve};
    const int j_c[3] = {pmb->jvs, pmb->jvs+pmb->block_size.nx2/2, pmb->jve};
    const int i_c[3] = {pmb->ivs, pmb->ivs+pmb->block_size.nx1/2, pmb->ive};

    // spatial cube nodes (k,j,i)~3x3x3 related to array indices x3
    const int cnr_idx[3][3][3][3] = {
      {{{kcs, jcs ,ics}, {kcs, jcs, icm}, {kcs, jcs, ice}},  // 0,0,0,:
       {{kcs, jcm, ics}, {kcs, jcm, icm}, {kcs, jcm, ice}},
       {{kcs, jce, ics}, {kcs, jce ,icm}, {kcs, jce, ice}}},
      {{{kcm, jcs ,ics}, {kcm, jcs, icm}, {kcm, jcs, ice}},  // 1,0,0,:

       {{kcm, jcm, ics}, {kcm, jcm, icm}, {kcm, jcm, ice}},
       {{kcm, jce, ics}, {kcm, jce ,icm}, {kcm, jce, ice}}},
      {{{kce, jcs ,ics}, {kce, jcs, icm}, {kce, jcs, ice}},  // 2,0,0,:
       {{kce, jcm, ics}, {kce, jcm, icm}, {kce, jcm, ice}},
       {{kce, jce, ics}, {kce, jce ,icm}, {kce, jce, ice}}}
    };

    // int cnr_idx[3][3][3][3];
    // for (int k=0; k<=2; ++k)
    //   for (int j=0; j<=2; ++j)
    //     for (int i=0; i<=2; ++i) {
    //       cnr_idx[k][j][i][2] = k_c[k];
    //       cnr_idx[k][j][i][1] = j_c[j];
    //       cnr_idx[k][j][i][0] = i_c[i];
    //     }

    // BD: debug [temp.] populate node_multiplicities
    for (int k=0; k<=2; ++k)
      for (int j=0; j<=2; ++j)
        for (int i=0; i<=2; ++i) {
          int ix_k = cnr_idx[k][j][i][0];
          int ix_j = cnr_idx[k][j][i][1];
          int ix_i = cnr_idx[k][j][i][2];
          // assumed +ve during dbh
          node_mult[k][j][i] = (int) (var(0, ix_k, ix_j, ix_i) + 0.5);
        }

    coutBoldRed("node_mult\n");
    print_int_arr3(&node_mult[0][0][0]);

    printf("AAA%d\n", *nblevel[0][0][0]);

    coutBoldRed("pmb->pbval->nblevel\n");
    print_int_arr3(&(pmb->pbval->nblevel[0][0][0]));

    int &lev = pmb->loc.level;
    coutBoldRed("lev =");
    printf("%d\n", lev);

    // apply corner conditions [top]
    //
    // N-W
    /*
    if (false)
    if (node_mult[0][0][0] > 1)
      if (node_mult[0][0][0] == 8) {
        // need to iterate in x,y,z
        for (int k=k_c[0]-ng; k<k_c[0]; ++k) {
          for (int ix=1; ix<ng; ++ix) {
            var(0, k, j_c[0]-ix, i_c[0]) /= 2.;
            var(0, k, j_c[0], i_c[0]-ix) /= 2.;
          }
          var(0, k, j_c[0], i_c[0]) /= 4.;
        }
        var(0, k_c[0], j_c[0], i_c[0]) /= 8.;
      } else if (node_mult[0][0][0] == 2) {
        // only node contributes
        var(0, ng, ng, ng) /= 2.;
      }
    */

    // Note:
    //
    // Various loop conditions can be written more compactly in what follows.
    // This would require multiple opers on entries which could give a race
    // cond. so we avoid this

    //-------------------
    // corner conditions
    //
    // 8,2,1 are unambiguous
    // 4 requires neighbour inspection
    // other multiplicities should be thrown as err
    //
    // N-W
    if (node_mult[0][0][0] != 1) {
      if (node_mult[0][0][0] == 8) {
        for (int ix=1; ix<=ng; ++ix) {
          var(0, k_c[0], j_c[0], i_c[0] - ix) /= 4.;
          var(0, k_c[0], j_c[0] - ix, i_c[0]) /= 4.;
          var(0, k_c[0] - ix, j_c[0], i_c[0]) /= 4.;
          for (int jx=1; jx<=ng; ++jx) {
            var(0, k_c[0], j_c[0] - jx, i_c[0] - ix) /= 2.;
            var(0, k_c[0] - jx, j_c[0] - ix, i_c[0]) /= 2.;
            var(0, k_c[0] - jx, j_c[0], i_c[0] - ix) /= 2.;
          }
        }
        var(0, k_c[0], j_c[0], i_c[0]) /= 8.;
      } else if (node_mult[0][0][0] == 2) {
        var(0, k_c[0], j_c[0], i_c[0]) /= 2.;
      } else if (node_mult[0][0][0] == 4) {

        if (pmb->pbval->nblevel[0][0][0] >= lev) {
          for (int jx=1; jx<=ng; ++jx) {
            for (int ix=1; ix<=ng; ++ix) {
              var(0, k_c[0], j_c[0]-jx, i_c[0]-ix) /= 2.;
            }
          }

          for (int ix=1; ix<=ng; ++ix) {
            var(0, k_c[0]-ix, j_c[0], i_c[0]) /= 2.;
          }

          var(0, k_c[0], j_c[0], i_c[0]) /= 4.;
        } else {

          if ((pmb->pbval->nblevel[0][0][1] >= lev) ||
              (pmb->pbval->nblevel[0][1][0] >= lev))
            for (int ix=1; ix<=ng; ++ix)
              var(0, k_c[0]-ix, j_c[0], i_c[0]) /= 2.;

          if ((pmb->pbval->nblevel[0][0][1] >= lev) ||
              (pmb->pbval->nblevel[1][0][1] >= lev))
            for (int ix=1; ix<=ng; ++ix)
              var(0, k_c[0], j_c[0]-ix, i_c[0]) /= 2.;

          if ((pmb->pbval->nblevel[0][1][0] >= lev) ||
              (pmb->pbval->nblevel[1][1][0] >= lev))
            for (int ix=1; ix<=ng; ++ix)
              var(0, k_c[0], j_c[0], i_c[0]-ix) /= 2.;

          var(0, k_c[0], j_c[0], i_c[0]) /= 4.;
        }


      } else {
        ErrorUnknownMultiplicity();
      }

    }


    //-------------------

    //-------------------
    // now treat as block without refinement with periodic BC
    /*
    if (false) {
    for (int k=pmb->kms; k<=pmb->kpe; ++k) {
      for (int i=pmb->ims; i<=pmb->ipe; ++i) {
        var(0, k, pmb->jvs, i) /= 2.;
        var(0, k, pmb->jve, i) /= 2.;
      }

      for (int j=pmb->jms; j<=pmb->jpe; ++j) {
        var(0, k, j, pmb->ivs) /= 2.;
        var(0, k, j, pmb->ive) /= 2.;
      }
    }

    for (int j=pmb->jms; j<=pmb->jpe; ++j) {
      for (int i=pmb->ims; i<=pmb->ipe; ++i) {
        var(0, pmb->kvs, j, i) /= 2.;
        var(0, pmb->kve, j, i) /= 2.;
      }
    }
    }
    */
    //-------------------

    // inspect only connectivity
    if (pmb->gid == -1) {
    // if ((node_mult[0][0][0] != 8) &&
    //     (node_mult[0][0][0] != 4) &&
    //     (node_mult[0][0][0] != 1) &&
    //     (node_mult[0][0][0] != 1)) {
      coutBoldRed("nm MB::UWIL gid = ");
      printf("%d\n", pmb->gid);
      Q();
    }//  else {
    //   return;
    // }

    if ((node_mult[0][0][0] == 22)
        // (node_mult[0][0][1] == 1) &&
        // (node_mult[0][1][1] == 1) &&
        // (node_mult[0][1][0] == 1)
        ) {
      var.print_all("%1.1f");
      coutBoldRed("nm MB::UWIL gid = ");
      printf("%d\n", pmb->gid);
      Q();
    }//  else {
    //   return;
    // }

    //-----------------
    // if (node_mult[0][0][0] != 4)
    // for (int k=0; k<=ng; ++k)
    //   for (int j=0; j<=pmb->jvs; ++j)
    //     for (int i=0; i<=pmb->ivs; ++i)
    //       if (!((std::abs(var(0, k, j, i) - 1.) < 0.001)
    //             || (std::abs(var(0, k, j, i)) < 0.001))) {
    //         var.print_all("%1.1f");
    //         coutBoldRed("\n\nerror found:");
    //         coutBoldRed("MB::UWIL gid = ");
    //         printf("%d\n", pmb->gid);
    //         Q();
    //       }
    // for (int k=pmb->kis; k<=pmb->kie; ++k)
    for (int k=0; k<=NGHOST; ++k)
      for (int j=0; j<=pmb->jvs; ++j)
        for (int i=0; i<=pmb->ivs; ++i)
          if (!((std::abs(var(0, k, j, i) - 1.) < 0.001)
                || (std::abs(var(0, k, j, i)) < 0.001))) {
            var.print_all("%1.1f");
            coutBoldRed("\n\nerror found:");
            coutBoldRed("MB::UWIL gid = ");
            printf("%d\n", pmb->gid);
            Q();
          }

    // _FinalizeVert3();
    // _FinalizeVertexConsistency3(0, 0, -1);
    // _FinalizeVertexConsistency3(0, 0, 1);

    // _FinalizeVertexConsistency3(0, -1, 0);
    // _FinalizeVertexConsistency3(0, 1, 0);

    // _FinalizeVertexConsistency3(-1, 0, 0);
    // _FinalizeVertexConsistency3(1, 0, 0);


    // int ng = NGHOST;
    // int b_nx1 = pmb->block_size.nx1;
    // int b_nx2 = pmb->block_size.nx2;
    // int b_nx3 = pmb->block_size.nx3;

    // int &lev = pmb->loc.level;

    // overall approach:
    // initially assume no refinement
    // correct afterwards if false


    // individual logic for faces hard-coded
    //
    // only one oxi is non-zero, thus with (e.g.)
    // ox3 = ox2 = 0, ox1 = 1
    // we can inspect the face-connected neighbor level with:
    // pmb->pbval->nblevel[ox3 + 1][ox2 + 1][ox1 + 1]

    // neighbor placement and level
    // int ox3 = -1, ox2 = 0, ox1 = 0;
    // int nb_lev = pmb->pbval->nblevel[ox3 + 1][ox2 + 1][ox1 + 1];

    // int si = pmb->iis, ei = pmb->iie;
    // int sj = pmb->jis, ej = pmb->jie;
    // int sk = pmb->kis, ek = pmb->kie;

    // if (ox1 < 0)
    //   si = ei = pmb->ivs;
    // else if (ox1 > 0)
    //   si = ei = pmb->ive;

    // if (ox2 < 0)
    //   sj = ej = pmb->jvs;
    // else if (ox2 > 0)
    //   sj = ej = pmb->jve;

    // if (ox3 < 0)
    //   sk = ek = pmb->kvs;
    // else if (ox3 > 0)
    //   sk = ek = pmb->kve;

    // int sgn_ox1 = (ox1 < 0) ? -1 : 1;
    // int sgn_ox2 = (ox2 < 0) ? -1 : 1;
    // int sgn_ox3 = (ox3 < 0) ? -1 : 1;

    // int uox3 = std::abs(ox3);
    // int uox2 = std::abs(ox2);
    // int uox1 = std::abs(ox1);

    // // infer internal corner idx
    // int I_cnr = ng + (ox1 + 1) * (b_nx1 / 2);
    // int J_cnr = ng + (ox2 + 1) * (b_nx2 / 2);
    // int K_cnr = ng + (ox3 + 1) * (b_nx3 / 2);

    // // reduce to 2d
    // int I = b_nx1 * uox1;
    // int J = b_nx2 * uox2;
    // int K = b_nx3 * uox3;

    // printf("(I,J,K) = (%d,%d,%d)\n", I, J, K);
    // printf("(I_cnr,J_cnr,K_cnr) = (%d,%d,%d)\n",
    //        I_cnr, J_cnr, K_cnr);
    //Q();

    // int N = 0;
    // for (int n_=nl_; n_<=nu_; ++n_)
    //   for (int k=sk; k<=ek; ++k)
    //     for (int j=sj; j<=ej; ++j)
    //       for (int i=si; i<=ei; ++i) {
    //         var(n_, k, j, i) /= 2.;
    //         // var(n_, k)
    //       }

  }


  if (false)
    if (pmb->block_size.nx3 > 1) {
      int ng = NGHOST;
      int bh_nx1 = pmb->block_size.nx1 / 2;
      int bh_nx2 = pmb->block_size.nx2 / 2;
      int bh_nx3 = pmb->block_size.nx3 / 2;

      for (int n=0; n<pbval_->nneighbor; n++) {
        // BD: debug
        coutMagenta("iterating over neighbors; (n, pbval_->nneighbor)=");
        printf("(%d, %d)\n", n, pbval_->nneighbor);
        //--

        NeighborBlock& nb = pbval_->neighbor[n];
        int ix_ox1 = nb.ni.ox1 + 1;
        int ix_ox2 = nb.ni.ox2 + 1;
        int ix_ox3 = nb.ni.ox3 + 1;

        // infer information about shared vertex multiplicity

        // common to same and finer
        if (mylevel <= nb.snb.level) {
          // neighbor at the same or finer level

          // retain the type of neighbor connect
          NeighborConnect nb_t = nb.ni.type;
          nb_type[ix_ox3][ix_ox2][ix_ox1] = nb_t;

          // retain information refinement
          nb_fi1[ix_ox3][ix_ox2][ix_ox1] = nb.ni.fi1;
          nb_fi2[ix_ox3][ix_ox2][ix_ox1] = nb.ni.fi2;

          // central, shared node of neighbor [edge or face]
          node_mult[ix_ox3][ix_ox2][ix_ox1] += 1;

          if (nb_t == NeighborConnect::face) {
            // two of {ox1, ox2, ox3} == 0
            // remainder non-zero
            printf("f: (ox1, ox2, ox3, l)=(%d, %d, %d, %d)\n",
                   nb.ni.ox1, nb.ni.ox2, nb.ni.ox3,
                   pmb->pbval->nblevel[ix_ox3][ix_ox2][ix_ox1]);

          } else if (nb_t == NeighborConnect::edge) {
            // one of ox1 == 0 | ox2 == 0 | ox3 == 0
            // remainder non-zero

            printf("e: (ox1, ox2, ox3, l)=(%d, %d, %d, %d)\n",
                   nb.ni.ox1, nb.ni.ox2, nb.ni.ox3,
                   pmb->pbval->nblevel[ix_ox3][ix_ox2][ix_ox1]);

            node_mult[ix_ox3][ix_ox2][ix_ox1] += 1;
            if (nb.ni.ox1 == 0) {
              node_mult[ix_ox3][ix_ox2][ix_ox1-1] += 1;
              node_mult[ix_ox3][ix_ox2][ix_ox1+1] += 1;
            }

          } else if (nb_t == NeighborConnect::corner) {
            //
            printf("c: (ox1, ox2, ox3, l)=(%d, %d, %d, %d)\n",
                   nb.ni.ox1, nb.ni.ox2, nb.ni.ox3,
                   pmb->pbval->nblevel[ix_ox3][ix_ox2][ix_ox1]);

            // node_mult[ix_ox3][ix_ox2][ix_ox1] += 1;
          }

        }

      }

      // coutBoldRed("MB::UWIL gid = ");
      // printf("%d\n", pmb->gid);
      // if(pmb->gid == 4)
      //   Q();

      // apply conditions
      for (int k=0; k<=2; k++) {
        for (int j=0; j<=2; j++) {
          for (int i=0; i<=2; i++) {
            int nblev = pmb->pbval->nblevel[k][j][i];
            //int nb_t = pmb->pbval->ni

            int ix_cnr = ng + i * bh_nx1;
            int jx_cnr = ng + j * bh_nx2;
            int kx_cnr = ng + k * bh_nx3;

            if (mylevel <= nblev) {
              int ox1 = i-1, ox2 = j-1, ox3 = k-1;
              // int sgn_ox1 = (ox1 < 0) ? -1 : 1;
              // int sgn_ox2 = (ox2 < 0) ? -1 : 1;
              // int sgn_ox3 = (ox3 < 0) ? -1 : 1;


            }

          }
        }
      }


      // int b_nx1 = pmb->block_size.nx1;
      // int b_nx2 = pmb->block_size.nx2;
      // int b_nx3 = pmb->block_size.nx3;

      // iterate over neighbors applying consistency conditions
      for (int n=0; n<pbval_->nneighbor; n++) {
        NeighborBlock& nb = pbval_->neighbor[n];
        int ox1 = nb.ni.ox1, ox2 = nb.ni.ox2, ox3 = nb.ni.ox3;

        // neighbor level and connect type
        // int (*ptr_nblev)[3][3][3] = &(pmb->pbval->nblevel);
        // int nblev = *ptr_nblev[ox3 + 1][ox2 + 1][ox1 + 1];
        int nblev = pmb->pbval->nblevel[ox3 + 1][ox2 + 1][ox1 + 1];

        NeighborConnect nb_t = nb.ni.type;

        // ignore coarser neighbors
        if (mylevel > nblev)
          continue;

        int ix_cnr = ng + (ox1 + 1) * bh_nx1;
        int jx_cnr = ng + (ox2 + 1) * bh_nx2;
        int kx_cnr = ng + (ox3 + 1) * bh_nx3;

        int si = pmb->iis, ei = pmb->iie;
        int sj = pmb->jis, ej = pmb->jie;
        int sk = pmb->kis, ek = pmb->kie;

        if (nb_t == NeighborConnect::face) {
          // two of {ox1, ox2, ox3} must be zero
          // categorize based on non-zero index

          if (ox1 != 0)
            si = ei = (ox1 < 0) ? pmb->ivs : pmb->ive;

          if (ox2 != 0)
            sj = ej = (ox2 < 0) ? pmb->jvs : pmb->jve;

          if (ox3 != 0)
            sk = ek = (ox3 < 0) ? pmb->kvs : pmb->kve;


        } else if (nb_t == NeighborConnect::edge) {
          // one of {ox1, ox2, ox3} must be zero
          if (ox1 == 0) {

            //.
          }
        } else {
          // .
        }
        if (mylevel <= nblev)
          for (int n_=nl_; n_<=nu_; ++n_) {
            for (int k=sk; k<=ek; ++k) {
              for (int j=sj; j<=ej; ++j) {
                for (int i=si; i<=ei; ++i) {
                  //var(n_, k, j, i) /= 2.;
                }
              }
            }
          }

        if (nb_t == NeighborConnect::edge) {
          // one of {ox1, ox2, ox3} must be zero

        }

        if (nb_t == NeighborConnect::corner) {
          // none of {ox1, ox2, ox3} can be zero
          int sgn_ox1 = (ox1 < 0) ? -1 : 1;
          int sgn_ox2 = (ox2 < 0) ? -1 : 1;
          int sgn_ox3 = (ox3 < 0) ? -1 : 1;


          for (int n_=nl_; n_<=nu_; ++n_) {
            for (int jx=1; jx<=ng; jx++) {
              for (int ix=1; ix<=ng; ix++) {
                // var(n_, kx_cnr, jx_cnr + sgn_ox2 * jx, ix_cnr + sgn_ox1 * ix) /= 2.;
                // var(n_, kx_cnr + sgn_ox3 * jx, jx_cnr + sgn_ox2 * ix, ix_cnr) /= 2.;
                // var(n_, kx_cnr + sgn_ox3 * jx, jx_cnr, ix_cnr + sgn_ox1 * ix) /= 2.;
              }
            }
          }


        }


      }

    } else if (pmb->block_size.nx2 > 1) {

      int ng = NGHOST;
      int bh_nx1 = pmb->block_size.nx1 / 2;
      int bh_nx2 = pmb->block_size.nx2 / 2;

      for (int n=0; n<pbval_->nneighbor; n++) {
        // BD: debug
        coutMagenta("iterating over neighbors; (n, pbval_->nneighbor)=");
        printf("(%d, %d)\n", n, pbval_->nneighbor);
        //--

        NeighborBlock& nb = pbval_->neighbor[n];
        int ox1 = nb.ni.ox1, ox2 = nb.ni.ox2;

        // infer information about shared vertex multiplicity

        // common to same and finer
        if (mylevel <= nb.snb.level) {
          // neighbor at the same or finer level

          // retain the type of neighbor connect
          nb_type[0][ox2 + 1][ox1 + 1] = nb.ni.type;

          // retain refinement information
          nb_fi1[0][ox2 + 1][ox1 + 1] = nb.ni.fi1;

          // central, shared node of neighbor [edge or face]
          node_mult[0][ox2 + 1][ox1 + 1] += 1;

          // north-south contributions to 2d corner
          if (ox1 == 0) {
            node_mult[0][ox2 + 1][0] += 1;
            node_mult[0][ox2 + 1][2] += 1;

            // finer level -> only one "sub-block"
            if (mylevel < nb.snb.level) {
              if (nb.ni.fi1 == 1) {
                node_mult[0][ox2 + 1][0] -= 1;
              } else {
                node_mult[0][ox2 + 1][2] -= 1;
              }
            }
          }

          // east-west
          if (ox2 == 0) {
            node_mult[0][0][ox1 + 1] += 1;
            node_mult[0][2][ox1 + 1] += 1;

            // finer level -> only one "sub-block"
            if (mylevel < nb.snb.level) {
              if (nb.ni.fi1 == 1) {
                node_mult[0][0][ox1 + 1] -= 1;
              } else {
                node_mult[0][2][ox1 + 1] -= 1;
              }
            }
          }

        }
      }

      // apply conditions

      for (int j=0; j<=2; j++) {
        for (int i=0; i<=2; i++) {
          int nblev = pmb->pbval->nblevel[1][j][i];

          int ix_cnr = ng + i * bh_nx1;
          int jx_cnr = ng + j * bh_nx2;

          // apply conditions on internal shared centers
          if (node_mult[0][j][i] > 1) {
            for (int n_=nl_; n_<=nu_; ++n_) {
              var(n_, 0, jx_cnr, ix_cnr) /= node_mult[0][j][i];
            }
          }

          if (mylevel <= nblev) {
            int ox1 = i-1, ox2 = j-1;
            int sgn_ox1 = (ox1 < 0) ? -1 : 1;
            int sgn_ox2 = (ox2 < 0) ? -1 : 1;

            if (nb_type[0][j][i] == NeighborConnect::edge) {
              for (int n_=nl_; n_<=nu_; ++n_) {
                for (int ix=1; ix<=ng; ix++) {
                  var(n_, 0, jx_cnr, ix_cnr + sgn_ox1 * ix) /= 2.;
                  var(n_, 0, jx_cnr + sgn_ox2 * ix, ix_cnr) /= 2.;
                }
              }
            } else if (nb_type[0][j][i] == NeighborConnect::face) {

              if (ox1 == 0) {
                for (int n_=nl_; n_<=nu_; ++n_) {
                  for (int ix=1; ix<bh_nx1; ix++) {
                    var(n_, 0, jx_cnr, ix_cnr - ix) /= 2.;
                    var(n_, 0, jx_cnr, ix_cnr + ix) /= 2.;
                  }
                }

              } else { // ox1 != 0 and ox2 = 0
                for (int n_=nl_; n_<=nu_; ++n_) {
                  for (int ix=1; ix<bh_nx2; ix++) {
                    var(n_, 0, jx_cnr - ix, ix_cnr) /= 2.;
                    var(n_, 0, jx_cnr + ix, ix_cnr) /= 2.;
                  }
                }

              }

              if (mylevel < nblev) {
                if (ox1 == 0) {
                  for (int n_=nl_; n_<=nu_; ++n_) {
                    for (int ix=1; ix<=ng; ix++) {
                      var(n_, 0, jx_cnr + sgn_ox2 * ix, ix_cnr) /= 2.;
                    }
                  }

                } else { // ox1 != 0 and ox2 = 0
                  for (int n_=nl_; n_<=nu_; ++n_) {
                    for (int ix=1; ix<=ng; ix++) {
                      var(n_, 0, jx_cnr, ix_cnr + sgn_ox1 * ix) /= 2.;
                    }
                  }

                }


              }

            }

          }

        }
      }


    } else {
      // 1d
      for (int n=0; n<pbval_->nneighbor; n++) {
        // Same and finer level neighbours additively unpack;
        // this leads to nodes with multiplicity

        NeighborBlock& nb = pbval_->neighbor[n];
        if (mylevel <= nb.snb.level) {
          if (nb.ni.ox1 < 0) {   // nb is leftward (correct on each cpt.)
            for (int n_=nl_; n_<=nu_; ++n_) {
              var(n_, 0, 0, pmb->ivs) /= 2.;
            }
          } else {               // nb is rightward
            for (int n_=nl_; n_<=nu_; ++n_) {
              var(n_, 0, 0, pmb->ive) /= 2.;
            }
          }

        }

      }

    }



  //
