'''
Created on May 29, 2020

@author: Vlad
'''

import os
import shutil

from magus_align.alignment_context import AlignmentContext
from magus_align.decompose.decomposer import decomposeSequences
from magus_align.merge.merger import mergeSubalignments
from magus_tools import external_tools
from magus_configuration import Configs
from magus_helpers import sequenceutils
from magus_tasks import task

'''
Alignments are treated as "tasks", units of work that are written out to task files and 
processed as threads and/or compute nodes become available. 
MAGUS tasks will recursively generate MAGUS tasks over large subsets and MAFFT tasks over smaller subsets.
'''

def mainAlignmentTask():    
    args = {"workingDir" : Configs.workingDir, "outputFile" : Configs.outputPath,
            "subalignmentPaths" : Configs.subalignmentPaths, "sequencesPath" : Configs.sequencesPath,
            "backbonePaths" : Configs.backbonePaths, "guideTree" : Configs.guideTree,
            "inputConstraint": Configs.inputConstraint}
    ######## added potential input constraint path to alignment context
    task = createAlignmentTask(args)
    task.submitTask()
    task.awaitTask()
    
def createAlignmentTask(args):
    return task.Task(taskType = "runAlignmentTask", outputFile = args["outputFile"], taskArgs = args)

# create a task to update the subalignment
def createUpdateSubalignmentTask(args):
    return task.Task(taskType = "runUpdateSubalignmentTask",
                     outputFile = args["outputFile"], taskArgs = args)

# run the "update subalignment" task
def runUpdateSubalignmentTask(**kwargs):
    aln = sequenceutils.readFromFastaOrdered(subalignmentPath)
    
    # remove taxa that are in constraintTaxa
    for taxon in constraintTaxa:
        if taxon in aln:
            aln.pop(taxon)
    
    # delete all-gap columns for aln and write to outputFile
    sequenceutils.cleanGapColumnsFromAlignment(aln, outputFile)

def runAlignmentTask(**kwargs):
    '''
    The standard MAGUS task: 
    decompose the data into subsets, align each subset, and merge the subalignments.
    '''
    
    with AlignmentContext(**kwargs) as context:
        if context.sequencesPath is not None:
            Configs.log("Aligning sequences {}".format(context.sequencesPath))
        
        decomposeSequences(context)
        if Configs.onlyGuideTree:
            Configs.log("Outputting only the guide tree, as requested..")
            shutil.copyfile(os.path.join(context.workingDir, "decomposition", "initial_tree", "initial_tree.tre"), context.outputFile)
            return
        
        alignSubsets(context)
        mergeSubalignments(context)

def alignSubsets(context):
    if len(context.subalignmentPaths) > 0:
        Configs.log("Subalignment paths already provided, skipping subalignments..")
        return
    
    Configs.log("Building {} subalignments..".format(len(context.subsetPaths)))
    subalignDir = os.path.join(context.workingDir, "subalignments")
    if not os.path.exists(subalignDir):
        os.makedirs(subalignDir)
        
    mafftThreshold = max(Configs.mafftSize, Configs.decompositionMaxSubsetSize, Configs.recurseThreshold)
    
    for file in context.subsetPaths:
        subset = sequenceutils.readFromFasta(file)
        subalignmentPath = os.path.join(subalignDir, "subalignment_{}".format(os.path.basename(file)))
        context.subalignmentPaths.append(subalignmentPath)
        
        if os.path.exists(subalignmentPath):
            Configs.log("Existing subalignment file detected: {}".format(subalignmentPath))       
             
        elif len(subset) <= mafftThreshold or not Configs.recurse:
            Configs.log("Subset has {}/{} sequences, aligning with MAFFT..".format(len(subset), mafftThreshold))            
            subalignmentTask = external_tools.buildMafftAlignment(file, subalignmentPath)
            context.subalignmentTasks.append(subalignmentTask)
            
        else:
            Configs.log("Subset has {}/{} sequences, recursively subaligning with MAGUS..".format(len(subset), mafftThreshold))
            subalignmentDir = os.path.join(subalignDir, os.path.splitext(os.path.basename(subalignmentPath))[0])
            subalignmentTask = createAlignmentTask({"outputFile" : subalignmentPath, "workingDir" : subalignmentDir, 
                                                    "sequencesPath" : file, "guideTree" : Configs.recurseGuideTree})   
            context.subalignmentTasks.append(subalignmentTask)

    task.submitTasks(context.subalignmentTasks)
    Configs.log("Prepared {} subset alignment tasks..".format(len(context.subalignmentTasks)))
    
    ######## for constrained MAGUS ########
    # before running mergeSubAlignments task, we need to update the subalignments
    # by removing any input constraint sequences, and making the input constraint
    # as a new subalignment.

    if context.inputConstraint:
        Configs.log("Detected an input constraint alignment: {}".format(
            context.inputConstraint) + ", using it as one of the subalignments.")
        
        # read taxon names from it
        constraint_unaln = sequenceutils.readFromFasta(context.sequencesPath,
                removeDashes=True)    
        context.constraintTaxa = [seq.tag for seq in constraint_unaln]
       
        # update current subalignments to remove any constraint taxa
        for subalignmentPath in context.subalignmentPaths:
            # overwriting existing subalignment path for now (avoiding
            # creating new files)
            updateSubalignmentTask = createUpdateSubalignmentTask(
                    {'subalignmentPath': subalignmentPath,
                     'outputFile': subalignmentPath,
                     'constraintTaxa': context.constraintTaxa})
            context.updateSubalignmentTasks.append(updateSubalignmentTask)
        task.submitTasks(context.updateSubalignmentTasks)
        # BOUNDARY CASE: removing taxa may have a chance of removing
        #                ALL taxa from a subalignment. For these subalignments
        #                remove them from the context list
        newSubalignmentPaths = []
        for subalignmentPath in context.subalignmentPaths:
            if os.stat(subalignmentPath).st_size > 0:
                newSubalignmentPaths.append(subalignmentPath)
        Configs.log("{}/{} original subalignments are modified and preserved".format(
            len(newSubalignmentPaths), len(context.subalignmentPaths)) + 
            " after removing input constraint taxa.")
        context.subalignmentPaths = newSubalignmentPaths

        # finally, add input constraint as one of the subalignment
        context.subalignmentPaths.append(context.inputConstraint)
        Configs.log("Added the input constraint alignment as one of the subalignment.")

