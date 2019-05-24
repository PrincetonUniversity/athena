By default, Athena++ assumes that every MeshBlock incurs roughly the same computational cost; therefore, the executable tries to distribute MeshBlocks as evenly as possible at the beginning of a simulation and when adaptive mesh refinement generates/destroys MeshBlocks. However, certain solver configurations involve additional physics that may cause uneven computational cost. Athena++ provides two ways to adjust the load balance other than the default one. 

*Note: these features may cause lower performance if they are not used correctly and there are some limitations, so please read the following description carefully.*

## Manual load balancing
If the cost of each MeshBlock is predictable and stable, one can use `MeshBlock::SetCostForLoadBalancing` in `MeshBlock::UserWorkInLoop` to assign the cost manually. Athena++ tries to redistribute MeshBlocks so that the total cost per node becomes as even as possible. This function must be called from all the MeshBlocks.
```c++
void MeshBlock::UserWorkInLoop() {
  Real mycost;
  // calculate the cost for this MeshBlock here
  SetCostForLoadBalancing(mycost);
  return;
}
```

And the manual load balancing method must be specified in the input file as follows:
```
    <loadbalancing>
    balancer   = manual    # load balancing method (default = "default")
    interval   = 10        # interval between load balancing (default = 10)
    tolerance  = 0.5       # acceptable load imbalance (default = 0.5 = 50%)
```
The interval parameter sets how frequently (in number of cycles) the load balancing is performed. Since the load balancing itself is somewhat costly as it involves global communication, we recommend performing it only occasionally. Note that when AMR is in use and one or more MeshBlocks are refined/derefined within a cycle, load balancing is immediately performed regardless of the interval parameter.

The tolerance parameter specifies the maximum permissible load imbalance. When the total cost in a node exceeds the average cost by this factor, it triggers the load balancing process. Because of the limited flexibility of the load balancing (see the notes below), this value should not be set too low. Also, for AMR, this parameter is ignored and set automatically because the adjustability changes according to the number of MeshBlocks (this implementation may subject to change).

## Automatic load balancing
The other option is automatic load balancing based on timing. With this feature enabled, the code measures the computation time for each MeshBlock and redistribute the load accordingly. To enable this feature, select the automatic load balancer in the input file:
```
    <loadbalancing>
    balancer   = automatic # load balancing method (default = "default")
    interval   = 10        # interval between load balancing (default = 10)
    tolerance  = 0.5       # acceptable load imbalance (default = 0.5 = 50%)
```
The other parameters remain the same as in the manual load balancing.

## Notes on load balancing
First of all, **the use of these load balancing features is recommended only when there is physics causing substantially uneven load**. The default load balancing is optimal when the cost is proportional to the number of MeshBlocks.

The load balancing works well only when there are sufficiently many MeshBlocks per node to allow flexible load redistribution. If there is one or two MeshBlocks per node, it is likely that the new load distribution does not improve the performance significantly, or it is even possible that no better solution can be found. Roughly speaking, the number of MeshBlocks per node needed for good load balancing should be as large as the cost ratio between the highest and lowest cost MeshBlocks. These features should be used only when the load imbalance is not negligible because smaller MeshBlocks can cause poor performance and you may achieve a better performance by ignoring minor load imbalance.

Also, it is not always possible to improve the performance because the code can redistribute MeshBlocks only along the Z-ordering. If the high cost MeshBlocks are badly localized on the Z-ordering, it is possible that the code cannot adjust the load sufficiently. This is the limitation of the current implementation.


## Adding a new task to automatic load balancing
When a new task is added, it should be registered to the load balancer if the cost of that task should be counted. This is controlled by `lb_flag` Boolean variable in the Task structure. If this flag is `true`, the code measures the time spent by this task, and it is ignored if the flag is `false`. It should be set in the `AddTask` member function of the `TaskList` class as follows:
```c++
    case (MY_PHYSICS):
      task_list_[ntasks].TaskFunc=
          static_cast<TaskStatus (TaskList::*)(MeshBlock*,int)>
          (&TimeIntegratorTaskList::MyPhysicsTask);
      task_list_[ntasks].lb_time = true;
      break;
```
Tasks involving MPI receiving should not be counted, as its time depends on the condition of the network and the behavior of other nodes.
