# omni_condition_timing.py
#
# Run the following
#
#   python omni_condition_timing.py -h
#
# for documentation, usage, and notes
#
# Takes a set of data tables each representing a part of a psychometric task
# Produces a set of .1D files containing trial timing information, one file for
# each type of trial, for each subject.
#
# By bmk27@cornell.edu, based on a script by Roy Moral, kms, and ebr

from pathlib import Path
import re
import csv
import ast

### Define some constants & defaults
TR_LENGTH = 2.5
TR_DROP = 4.0 # MAKE CERTAIN THIS IS SET THE SAME AS MEICA_BATCH_CLUSTER.SH
ZERO_TIME = TR_DROP*TR_LENGTH
EXPECTED_RUN_NUMBERS = [1, 2, 3, 4]
# Prefix & extension that input data table files will have
DATA_PREFIX = 'tempAttnAudTME'
DATA_EXTENSION = '.txt'
# Which post-task columns will be merged into the mid-task data?
POST_TASK_MERGE_COLUMNS = ['acc', 'oldLoc', 'newLoc', 'confRT', 'confVal', 'oldName', 'newName']
# Prefix and extension to use when creating output timing files
OUTPUT_PREFIX = 'stimes'
OUTPUT_EXTENSION = '.1D'

def find_data_files(mid_task_data_directory, post_task_data_directory, data_prefix='tempAttnAudTME', data_extension='.txt'):
    # find_data_files:
    #   Find all mid-task and post-task data files, and organize them by subject ID
    #   Arguments:
    #       mid_task_data_directory = path to the directory where the mid task data files can be found
    #       post_task_data_directory = path to the directory where the mid task data files can be found
    #       data_prefix = (optional) the string prefix for the data files
    #       data_extension = (optional) the string extension for data files
    #   Returns:
    #       A dictionary of data files and metadata, organized by subject ID

    # Ensure data extension begins with a dot
    if data_extension[0] != '.':
        data_extension = '.' + data_extension
    # Ensure directories are Path objects
    mid_task_data_directory = Path(mid_task_data_directory)
    post_task_data_directory = Path(post_task_data_directory)

    # Find files matching pattern
    mid_task_data_files = mid_task_data_directory.glob('{prefix}*_ef{extension}'.format(prefix=data_prefix, extension=data_extension))
    post_task_data_files = post_task_data_directory.glob('{prefix}*{extension}'.format(prefix=data_prefix, extension=data_extension))
    data_files = set(list(mid_task_data_files) + list(post_task_data_files))
    # Loop over data files, extract the subject ID for each, and
    subjects = {}
    for data_file in data_files:
        # Parse filename
        ID, age, gender, run, type = get_file_info(data_file, data_prefix=data_prefix)
        if ID is None:
            # Data file didn't match pattern. Skip it.
            continue
        # Initialize subject dictionary if this is the first file for them
        if ID not in subjects:
            subjects[ID] = {'ID':ID, 'age':age, 'gender':gender, 'run_files':{}, 'post_file':None, 'warnings':[], 'errors':[]}

        if type == 'mid_task':
            # This is a mid-task data file from one of the four runs
            # Make sure we're not finding multiple mid-task data files for the same subject & run
            if run in subjects[ID]['run_files']:
                subjects[ID]['errors'].append('Error, two different mid-task data files for subject #{ID}, run #{run} have been found. That should not happen.'.format(ID=ID, run=run))
            # Store mid-task data file under the correct subject and run
            subjects[ID]['run_files'][run] = data_file
        elif type == 'post_task':
            # This is a post-task assessment data file
            # Make sure we're not finding multiple post-task data files for the same subject
            if subjects[ID]['post_file'] is not None:
                subjects[ID]['errors'].append('Error, two different post-task data files for subject #{ID} have been found. That should not happen.'.format(ID=ID))
            # Store post-task data file under the correct subject
            subjects[ID]['post_file'] = data_file

    # Do some sanity checking to make sure we have a full, valid data set for each subject
    for ID in subjects:
        if len(subjects[ID]['run_files']) == 0:
            # Zero mid-task files found
            subjects[ID]['errors'].append('No mid-task files found for subject {ID}. Dropping subject.'.format(ID=ID))
        else:
            for run_num_expected in EXPECTED_RUN_NUMBERS:
                if run_num_expected not in subjects[ID]['run_files'].keys():
                    subjects[ID]['warnings'].append('Subject {ID} missing data file for run {k}.'.format(ID=ID, k=run_num_expected))
                    subjects[ID]['run_files'][run_num_expected] = None
        if subjects[ID]['post_file'] is None:
            subjects[ID]['errors'].append('No post-task file found for subject {ID}. Dropping subject.'.format(ID=ID))

    return subjects

def get_file_info(data_file, data_prefix='tempAttnAudTME'):
    # get_file_info:
    #   Parse subject and run info from data filename, for example tempAttnAudTME-14.f21-1.1.l.1_ef.txt
    #   Arguments:
    #       data_file = path to a data file
    #       data_prefix  = the string prefix for data files
    #   Returns:
    #       tuple of (ID, age, gender, run, type) = a list of information
    #           extracted from the data file name

    # Returns subject ID, subject gender, subject age, and run #
    data_file_name = Path(data_file).name
    pattern = '{prefix}\-([0-9]+)\.([mf])([0-9]+)\-([0-9]+)\.([0-9]+)'.format(prefix=data_prefix)
    m = re.match(pattern, data_file_name)
    if m is None:
        return None, None, None, None, None

    [ID, gender, age, type, run] = m.groups()
    ID = int(ID)
    age = int(age)
    run = int(run)
    type = int(type)
    # Interpret data file type
    if type == 1:
        type = 'mid_task'
    elif type == 2:
        type = 'post_task'
    return ID, age, gender, run, type

def load_data_file(data_file):
    # load_file_info:
    #   Read data file into a dictionary of columns
    #   Arguments:
    #       data_file = path to a data file
    #   Returns:
    #       data = a dictionary of columns representing the data in data_file

    with open(data_file) as f:
        # Read tab-delimited info
        reader = csv.reader(f, delimiter="\t")
        # Convert to list of row-lists
        raw_data = list(reader)
    # Make sure file is not empty
    if len(raw_data) == 0:
        return None
    # Extract header row
    header = raw_data.pop(0)
    # Reorganize into a dictionary where each key is a heading, and each value is a list of values for each trial for that heading
    data = {}
    for col, field in enumerate(header):
        data[field] = convert_column_type([row[col] for row in raw_data])
    return data

def add_run_field(data, run):
    # add_run_field:
    #   Add data column called "run" in which each entry is the run number
    #   Arguments:
    #       data = dictionary of data columns (keys are field names)
    #       run = integer run number
    #   Returns:
    #       data = a new dictionary of data columns including the new "run"
    #           field
    #
    # Get any field to find out how long the run column should be
    example_field = list(data.keys())[0]
    # Generate run column
    data['run'] = [run for datum in data[example_field]]
    return data

def merge_post_task_data(mid_task_data, post_task_data,
                    post_task_merge_columns=POST_TASK_MERGE_COLUMNS,
                    mid_task_merge_key='imgName', post_task_merge_key='oldName'):
    # merge_post_task_data:
    #   Merge specified columns from post-task data into mid-task data
    #   Arguments:
    #       mid_task_data = dictionary of mid-task data
    #       post_task_data = dictionary of post-task data
    #       post_task_merge_columns = (optional) a list of names of fields to
    #           merge from post-task data into mid-task data. Any other fields
    #           in post-task data will be ignored, and will not be available
    #           for use as a condition
    #       mid_task_merge_key = (optional) column in the mid-task data files
    #           that will be used to merge data
    #       post_task_merge_key = (optional) column in the post-task data files
    #           that will be used to merge data
    #   Returns:
    #       mid_task_data = dictionary of mid-task data with specified
    #           post-task data merged into it

    for post_task_column in post_task_merge_columns:
        if post_task_column in mid_task_data:
#            print('    Merging post-task {col} into mid-task data, but {col} already exists in mid-task data. Deleting {col} from mid-task data...'.format(col=post_task_column))
            del mid_task_data[post_task_column]

    # Create empty columns in mid-task data to hold merged data
    mid_task_data_len = len(mid_task_data['subid'])
    for new_column in post_task_merge_columns:
        mid_task_data[new_column] = [None for k in range(mid_task_data_len)]

    # Loop over merge column in mid task data
    mid_task_data_len = len(mid_task_data[mid_task_merge_key])
    for mid_task_row, merge_key in enumerate(mid_task_data[mid_task_merge_key]):
        # Loop over merge column in post task data and look for a match
        for post_task_row, potential_merge_key in enumerate(post_task_data[post_task_merge_key]):
            if merge_key == potential_merge_key:
                # We have a match! Merge the specified columns in
                for new_column in post_task_merge_columns:
                    mid_task_data[new_column][mid_task_row] = post_task_data[new_column][post_task_row]
                # Each merge key will only appear once in the post task data, so we can move on
                break

    # Return mid_task data with post_task data merged in
    return mid_task_data

def convert_column_type(data_column):
    # convert_column_type:
    #   Take a column of data values, and attempt to parse it as an
    #       integer or float
    #   Arguments:
    #       data_column = a list of data values
    #   Returns:
    #       data_column = a list of data values, converted to int or float if
    #           possible.

    # Loop over potential parsers to see which one successfully converts entire column.
    parsers = [int, float]
    for parser in parsers:
        new_data_column = []
        parser_worked = True
        # Loop over column, try parser on each element.
        for datum in data_column:
            try:
                # Try to parse datum
                new_data_column.append(parser(datum))
            except ValueError:
                # Parser did not work. Move on to next parser.
                parser_worked = False
                break
        if parser_worked:
            # Parser successfully converted entire column without error
            return new_data_column
    return data_column

def get_timestamps(subject_data, condition_fields, condition_classes, timing_field='cycleOnset', concatenate_runs=False, concatenation_timing_field='nextOnset'):
    # get_timestamps:
    #   Gather timestamps for given condition fields for all subjects
    #   Arguments:
    #       subject_data = a dictionary of info about subject data
    #       condition_fields = a list of names of fields to use as conditions
    #       condition_classes = a dictionary of classes to group values into,
    #           one for each field. If None is passed in for a field instead
    #           of a dictionary, then that field will not be grouped.
    #           For example, [{'type1':[0, 1, 2], 'type2':[3, 4, 5]}, ...] would
    #           rename any values equal to 0, 1, or 2 to "type1", and any values
    #           equal to 3, 4, or 5 to "type2". If this entire argument is None,
    #           no fields will be grouped.
    #       concatenate_runs = a boolean flag indicating that multiple runs
    #           should be combined into a single series of timestamps. Since
    #           timestamps in subsequent runs reset to zero, the
    #           concatenation_timing_field will be used to determine the
    #           starting timestamp of each run following run #1.
    #           If this option is True, subjects with runs that don't begin
    #           with #1 and continue consecutively will be dropped with fatal
    #           errors
    #       concatenation_timing_field = the field that will be used to join
    #           the timestamps of subsequent runs ONLY if concatenate_runs is
    #           True. Otherwise this argument is ignored.
    #   Returns:
    #       timestamps = a dictionary of timestamp data organized by subject ID
    #           then run number
    #       warnings = a list of warning messages, if any non-fatal errors occurred
    #       error = an error message, if a fatal error occurred

    # Load post-task data, for merging in relevant fields

    warnings = []

    post_task_data_file = subject_data['post_file']
    post_task_data = load_data_file(post_task_data_file)
    if post_task_data is None:
        # Something wrong with the data file, we'll skip this subject.
        return None, warnings, 'Post task data file for subject {ID} was empty, missing, or otherwise invalid.'.format(ID=subject_data['ID'])

    timestamps = {}
    count = 0
    runs = sorted(subject_data['run_files'].keys())
    run_end_times = [None for run_idx in EXPECTED_RUN_NUMBERS]
    skipped_run = False
    for run_idx, run in enumerate(runs):
        # Get mid-task data file path
        data_file = subject_data['run_files'][run]

        # Handle situation if run file is missing
        if data_file is None:
            # Missing run file
            skipped_run = True
            warnings.append('Subject {ID} missing mid-task file for run {n}.'.format(ID=subject_data['ID'], n=run))
            continue
        elif concatenate_runs:
            if skipped_run == True:
                # Problem: We have a skipped run, followed by an unskipped run,
                #   and concatenate_runs is True. We can't concatenate
                #   subsequence runs with missing runs. Drop the subject.
                return None, warnings, 'Subject {ID} has files for runs that are subsequent to missing runs, and run concatenation has been requested. Cannot concatenate runs if an earlier run is missing - dropping subject.'.format(ID=subject_data['ID'])

        # Load mid-task data
        data = load_data_file(data_file)

        # Handle situation if data file can't be loaded
        if data is None:
            # Something wrong with the data file, we'll skip this subject.
            warnings.append('Run {n} mid task data file for subject {ID} was empty, missing, or otherwise invalid.'.format(ID=subject_data['ID'], n=run))
            continue

        # Merge post-task data in
        data = merge_post_task_data(data, post_task_data)

        # Add synthetic "run" field into data so it's possible to select by run number
        data = add_run_field(data, run)

        # Ensure all condition fields exist as fields in the merged data table
        for condition_field in condition_fields:
            if condition_field not in data:
                raise KeyError('Condition {c} not found in data table. Please check spelling/capitalization and try again.'.format(c=condition_field))

        # If concatenate_runs is True, we'll need the end times for each run so we can add them on to subsequent runs
        if concatenate_runs:
            run_end_times[run_idx] = data[concatenation_timing_field][-1]

        # Apply value groupings, if any are requested
        data = apply_value_grouping(data, condition_fields, condition_classes)

        # Generate a list of all conditions for the requested condition_fields
        conditions = list(zip(*[data[condition_field] for condition_field in condition_fields]))
        # Generate a set of unique conditions for the requested condition_fields
        unique_conditions = set(conditions)
        start_timestamp = data[timing_field][0]
        for target_condition in unique_conditions:
            # Create timing file for this condition
            if target_condition not in timestamps:
                # Create blank lists to hold timestamps
                timestamps[target_condition] = [[] for run in EXPECTED_RUN_NUMBERS]
#            timestamps[target_condition].append([])
            foundCondition = False
            for condition, timestamp in zip(conditions, data[timing_field]):
                if condition == target_condition:
                    foundCondition = True
                    # Adjust timestamps
                    if concatenate_runs and run_idx > 0:
                        # We're concatenating runs, so we need to add on the cumulative ending times for previous runs
                        new_timestamp = timestamp - start_timestamp - ZERO_TIME + sum(run_end_times[:run_idx])
                    else:
                        # Not concatenating_runs, proceed as usual
                        new_timestamp = timestamp - start_timestamp - ZERO_TIME
                    # Record timestamp
                    timestamps[target_condition][run_idx].append(new_timestamp)
                    count = count + 1
            if not foundCondition:
                print('Found no match for ', target_condition)

    if concatenate_runs:
        # Concatenate run timestamps
        for target_condition in timestamps:
            concatenated_timestamps = []
            for run_timestamps in timestamps[target_condition]:
                if run_timestamps is not None:
                    concatenated_timestamps.extend(run_timestamps)
            timestamps[target_condition] = [concatenated_timestamps]

    return timestamps, warnings, None

def apply_value_grouping(data, condition_fields, condition_classes):
    # apply_value_grouping
    #   Arguments:
    #       data = dictionary containing a column of values for each field key
    #       condition_fields = a list of names of fields to use as conditions
    #       condition_classes = a dictionary of classes to group values into,
    #           one for each field. If None is passed in for a field instead
    #           of a dictionary, then that field will not be grouped.
    #           For example, [{'type1':[0, 1, 2], 'type2':[3, 4, 5]}, ...] would
    #           rename any values equal to 0, 1, or 2 to "type1", and any values
    #           equal to 3, 4, or 5 to "type2". If this entire argument is None,
    #           no fields will be grouped.
    if condition_classes is None:
        # No groupings to apply - return the data table unmodified
        return data

    # Loop over condition fields
    for condition_field, condition_class_spec in zip(condition_fields, condition_classes):
        if condition_class_spec is not None:
            # Loop over classes within grouping specification
            for class_name in condition_class_spec:
                # Loop over each value in the data column and replace it with its class, thus achieving grouping
                for k, value in enumerate(data[condition_field]):
                    if value in condition_class_spec[class_name]:
                        # Replace data value with the class name itself
                        data[condition_field][k] = class_name
    return data

def evaluate_condition_class_spec(condition_class_spec_string):
    # evaluate_condition_class_spec
    #   Convert a string condition class grouping specification into a python
    #       dictionary
    # Arguments:
    #       condition_class_spec_string = a user-supplied string representing a
    #           grouping scheme for the values for a data table field
    if condition_class_spec_string is None:
        return condition_class_spec_string

    format_good = True
    format_error = ''
    try:
        condition_class_spec = ast.literal_eval(condition_class_spec_string)
        if type(condition_class_spec) is not dict:
            # The grouping specification is not a dictionary
            format_good = False
            format_error = 'grouping specification is not a dictionary'
        elif not all([type(condition_class_spec[className]) in (list, tuple) for className in condition_class_spec]):
            # The groupings are not all python lists/tuples
            format_good = False
            format_error = 'groups are not all lists/tuples'
    except (SyntaxError, ValueError) as e:
        # Could not parse specification
        format_good = False
        format_error = 'grouping specification had invalid syntax'

    if not format_good:
        raise ValueError('Invalid condition class statement: {s} *** Error={e} *** Condition class statement must be a dictionary of the form {{class1:[val1, val2, ..., vallN], class2:[val1, val2, ..., valN], ..., classN:[val1, val2, ..., valN]}}'.format(s=condition_class_spec_string, e=format_error))

    return condition_class_spec

def create_condition_file(output_directory, subject_ID, condition,
    condition_fields, timestamps, prefix=OUTPUT_PREFIX,
    extension=OUTPUT_EXTENSION, concatenate_runs=False, skip_empty_runs=False):
    # create_condition_file:
    #   Using found timestamps, create a condition timing file containing the
    #       timestamps of any trials that match the given condition
    #   Arguments:
    #       output_directory = path to the directory in which to create the file
    #       subject_ID = the subject ID of the timing file (for constructing the filename)
    #       condition = the condition values for this timing file (for constructing the filename)
    #       condition_fields = the names of the condition fields (for constructing the filename)
    #       timestamps = a list of lists representing the timestamps for the
    #           given condition for each run for this subject
    #       prefix = (optional) the prefix to put at the beginning of the filename
    #       extension = (optional) the extension to append to the end of the filename
    #       concatenate_runs = a boolean flag indicating that multiple runs
    #           should be combined into a single series of timestamps. The suffix
    #           '-runconcat' will be added to the filename before the extension.
    #       skip_empty_runs = a boolean flag indicating whether or not runs with no times
    #           should get an empty row in the timing file. If True, the timing file will
    #           not have any blank rows. If False, the timing file will contain one or more
    #           blank rows if one or more runs have no times for the given condition
    #   Returns:
    #       timestamps = a dictionary of timestamp data organized by subject ID
    #           then run number

    # Ensure extension begins with a dot
    if extension[0] != '.':
        extension = '.' + extension

    # Construct the list of condition field/values for the filename
    condition_string = '-'.join([str(condition_field) + '_' + str(condition_value) for condition_field, condition_value in zip(condition_fields, condition)])

    if concatenate_runs:
        concat_string = '-runconcat'
    else:
        concat_string = ''

    # Construct the full filename
    file_name = '{prefix}-{ID}-{condition}{concat}{ext}'.format(prefix=prefix, ID=subject_ID, condition=condition_string, ext=extension, concat=concat_string)

    # Add on the path
    file_path = Path(output_directory) / Path(file_name)

    # Write the timestamp data to the file
    with open(file_path, 'w') as f:
        f.write('\n'.join([' '.join([str(t) for t in run_timestamps]) for run_timestamps in timestamps if ((not skip_empty_runs) or len(run_timestamps) > 0)]))

def print_help():
    # print_help:
    #   Print out usage help info
    print('omni_condition_timing:')
    print('    Create timing files across a given set of conditions, using both mid-task')
    print('    data and post-task data.')
    print()
    print('    Usage:')
    print()
    print('        python omni_condition_timing.py [mid-task-dir] [post-task-dir] field1')
    print('            field2 ... fieldN [output-dir]')
    print()
    print('            mid-task-dir:          Path to directory where mid-task data can be')
    print('                                   found')
    print('            post-task-dir:         Path to directory where post-task data can')
    print('                                   be found')
    print('            field1...N:            One or more field names to use as conditions.')
    print('                                   Each must match a column name in the')
    print('                                   mid_task or post_task data tables. Field')
    print('                                   names may be optionally followed by class ')
    print('                                   grouping statements. See note below.')
    print('                                   You may also use use the field name \'run\'')
    print('                                   to use run number as a condition.')
    print('            output-dir:            Path to directory where the output timing')
    print('                                   files should be written')
    print('            -c, --concatenate-runs Optional flag indicating that the')
    print('                                   timestamps for each run should not be reset')
    print('                                   at the beginning of each run, and should all')
    print('                                   be printed on a single line in the output')
    print('                                   file')
    print('            -s, --skip-empty-runs  Optional flag indicating that empty runs')
    print('                                   should NOT result in an empty row in the')
    print('                                   timing file. By default, a run that does not')
    print('                                   contain any times for the selected condition')
    print('                                   will produce an empty row. If you include this')
    print('                                   flag, empty rows will not be produced.')
    print()
    print('    Example:')
    print()
    print('        python omni_condition_timing.py ./mid-task ./post-task toneType ')
    print('            imgType ./timing')
    print()
    print('    Notes:')
    print('        ****Grouping classes****')
    print('        Each field name may be followed by a class grouping statement, ')
    print('        indicating that values of the field in question should be grouped ')
    print('        into classes, rather than treated as separate values.')
    print('        For example:')
    print('            python omni_condition_timing.py ./mid-task ./post-task toneType')
    print('                \'imgType:{"scrambled":[0],"unscrambled":[1,2,3,4,5,6,7]}\'')
    print('                ./timing')
    print('        Note that each class grouping statement is attached to its field name')
    print('        by a colon, and the statement that follows the colon is a valid')
    print('        python dictionary literal of the form:')
    print('        {class1:[val1,val2,...,vallN],class2:[val1,val2,...,valN],...,classN:[val1,val2,...,valN]}')
    print('        Note that the field argument with a class grouping statement needs')
    print('        to be enclosed in quotes that are not the same kind of quotes used')
    print('        within the class grouping statement, otherwise the shell may')
    print('        improperly parse the argument list.')
    print('        ****Run concatenation****')
    print('        Timestamps for Run N will be modified by adding on the sum of the ')
    print('        last values of the field "nextOnset" for runs 1 to N-1.')
    print('        Note that if run concatenation is on, and subjects have run N, but')
    print('        are missing run M such that M < N, then the subject will be dropped.')

if __name__ == '__main__':
    # This runs when the script is directly invoked (as opposed to imported as a library)
    import sys

    ### Parse command line arguments
    # Search for boolean flags first
    concatenate_runs = False
    skip_empty_runs = False
    flag_idx = []
    for idx, arg in enumerate(sys.argv):
        if arg == '-h':
            flag_idx.append(idx)
            print_help()
            exit()
        elif arg in ['-s', '--skip-empty-runs']:
            flag_idx.append(idx)
            skip_empty_runs = True
        elif arg in ['-c', '--concatenate-runs']:
            flag_idx.append(idx)
            concatenate_runs = True

    # Remove boolean flags from argument list
    sys.argv = [arg for idx, arg in enumerate(sys.argv) if idx not in flag_idx]

    num_args = len(sys.argv)
    arg_num = 0
    condition_fields = []
    condition_classes = []
    # Proceed to parse the rest of the arguments now that boolean flags have been removed.
    while len(sys.argv) > 0:
        arg = sys.argv.pop(0)
        if arg_num == 0:
            pass
        elif arg_num == 1:
            # Arg #1 is the directory where mid-task data files can be found
            mid_task_data_directory = Path(arg)
        elif arg_num == 2:
            # Arg #2 is the directory where post-task data files can be found
            post_task_data_directory = Path(arg)
        elif arg_num == num_args-1:
            # Last argument is output directory
            output_directory = Path(arg)
        elif arg_num > 2:
            # All remaining args are names of fields to be used as conditions
            if ':' in arg:
                # Grouping classes were provided for this field
                field, class_spec = arg.split(":", 1)
                class_spec = evaluate_condition_class_spec(class_spec)
            else:
                # No grouping classes provided
                field = arg
                class_spec = None
            condition_fields.append(field)
            condition_classes.append(class_spec)
        arg_num = arg_num + 1

    ### Sanity check inputs:
    if not mid_task_data_directory.is_dir():
        raise IOError('Error - given mid-task data directory is not a valid path to a directory: {path}. Run this script with the -h flag to get usage help.'.format(path=mid_task_data_directory))
    if not post_task_data_directory.is_dir():
        raise IOError('Error - given post-task data directory is not a valid path to a directory: {path}. Run this script with the -h flag to get usage help.'.format(path=post_task_data_directory))
    if not output_directory.exists():
        print('Given output directory not found - creating it:')
        print('    ' + str(output_directory))
        output_directory.mkdir()
    if len(condition_fields) == 0:
        raise ValueError('Please specify at least one condition field for which to find timing info. Run this script with the -h flag to get usage help.')

    ### Do it
    print()
    print('Finding subject files...')
    subjects = find_data_files(mid_task_data_directory, post_task_data_directory, data_prefix=DATA_PREFIX, data_extension=DATA_EXTENSION)
    print('...found {n} subjects'.format(n=len(subjects)))
    print()
    for k, subject_ID in enumerate(subjects):
        print('Processing subject {ID} ({pct}%)...'.format(ID=subject_ID, pct=int(100*k/len(subjects))))
        if len(subjects[subject_ID]['errors']) > 0:
            # Some kind of fatal error was found with this subject, don't continue processing.
            pass
        else:
            # Collect and organize timestamps for each subject, run, and condition
            timestamps, warnings, error = get_timestamps(subjects[subject_ID], condition_fields, condition_classes, concatenate_runs=concatenate_runs)
            if error is not None:
                # Something must be invalid - skip this subject
                subjects[subject_ID]['errors'].append(error)
            subjects[subject_ID]['warnings'].extend(warnings)
            if timestamps is not None:
                # Loop over detected conditions
                for condition in timestamps:
                    # Save condition timing file
                    create_condition_file(output_directory, subject_ID,
                        condition, condition_fields, timestamps[condition],
                        concatenate_runs=concatenate_runs,
                        skip_empty_runs=skip_empty_runs)
        # Print out warning & error messages for this subject
        for warning in subjects[subject_ID]['warnings']:
            print('    WARNING: ', warning)
        for error in subjects[subject_ID]['errors']:
            print('    ERROR:   ', error)

    # Classify subjects by which ones have warnings/errors
    perfectSubjects = [ID for ID in subjects if len(subjects[ID]['errors']) == 0 and len(subjects[ID]['warnings']) == 0]
    goodSubjects =    [ID for ID in subjects if len(subjects[ID]['errors']) == 0 and len(subjects[ID]['warnings']) > 0]
    badSubjects =     [ID for ID in subjects if len(subjects[ID]['errors']) > 0]
    print('...done processing subjects!')
    print()

    ### Print summary report
    print('******************** SUMMARY REPORT ********************')
    print('  Condition fields selected:')
    for condition_field, condition_class_spec in zip(condition_fields, condition_classes):
        if condition_class_spec is None:
            condition_class_text = ''
        else:
            condition_class_text = ' with grouping scheme: {c}'.format(c=condition_class_spec)
        print('    {c}'.format(c=(condition_field + condition_class_text)))
    print('  Run concatenation: {onoff}'.format(onoff=('on' if concatenate_runs else 'off')))
    print('  ')
    print('  Subjects without errors or warnings ({k} of {n}):'.format(k=len(perfectSubjects), n=len(subjects)))
    print('    ', ', '.join([str(ID) for ID in perfectSubjects]))
    print('  Subjects without errors, but with some warnings ({k} of {n}):'.format(k=len(goodSubjects), n=len(subjects)))
    print('    ', ', '.join([str(ID) for ID in goodSubjects]))
    print('  Subjects that failed to run because of errors ({k} of {n}):'.format(k=len(badSubjects), n=len(subjects)))
    print('    ', ', '.join([str(ID) for ID in badSubjects]))
