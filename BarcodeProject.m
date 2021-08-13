
%% BarcodeProject
FLAG_CALCULATE_DISTANCE_MATRIX = 1; % 0: load, 1: calculate

%% define problem
% possible nucleotides
Nuc = ['A','C','G','T'];
barcode_length = 4;

sequences = permn(Nuc,barcode_length);

set_size    = length(sequences); % number of possible sequences
subset_size = 8; % number of chosen sequences

%% calculate distance matrix
if FLAG_CALCULATE_DISTANCE_MATRIX == 1
distance_matrix = zeros(set_size,set_size);
for i = 1:set_size
    for j = 1:set_size
        distance_matrix(i,j) = strdist(sequences(i,:),sequences(j,:));
    end
end
max_dist = max(distance_matrix(:));
distance_matrix(distance_matrix==0) = max_dist+1;

%% load distance matrix
else
filename = 'kmer_7_distance_matrix.csv';
file = importdata(filename);
distance_matrix = file.data;

filename = 'kmer_6_distance_matrix.mat';
max_dist = max(distance_matrix(:));
distance_matrix(distance_matrix==0) = max_dist+1;
end

%% get random subset
subset_indices      = datasample(1:set_size,subset_size,'Replace',false)';
subset_indices_orig = subset_indices;

%% optimize subset
MAX_ITER = 200; ITERATIONS = 1;
max_iter = 100;
figure,hold on
for ITERATIONS = 1:MAX_ITER
    % find point with minimum distance within subset
    distance_matrix_subset = distance_matrix(subset_indices,subset_indices);
    max_dist = max(distance_matrix_subset(:));
    distance_matrix_subset(distance_matrix_subset==0) = max_dist + 1;
    
    % minimum distance from each point individually
    [D_old_vec,distance_matrix_index_old_vec] = min(distance_matrix_subset);
    
    % minimum and mean distances within subset
    [D_old,index_old] = min(D_old_vec);
    mean_D_old = mean(D_old_vec);
    
    % get indices of remaining points
    point_old = subset_indices(index_old);
    query_subset_indices = subset_indices;
    query_subset_indices(index_old) = [];
    
    D_new = D_old; % initialize distances
    for iterations = 1:max_iter
        % get new random point  
        set_to_sample_from = 1:set_size;
        for s = 1:subset_size
            set_to_sample_from =  set_to_sample_from(set_to_sample_from ~= subset_indices(s));
        end
        index_new = datasample(set_to_sample_from,1,'Replace',false)';
        
        % get distance of new point to query subset
        distance_vector = distance_matrix(index_new,query_subset_indices);
        mean_D_new = mean(distance_vector);
        D_new      = min(distance_vector);
        
        % exit loop if better point has been found
        if D_old < D_new
            break
        end
        
        if D_old == D_new
            mean_distance_old = (sum(sum(distance_matrix(subset_indices,subset_indices)))-(subset_size*max_dist))./(subset_size*subset_size-subset_size);
            mean_distance_new = (sum(sum(distance_matrix([query_subset_indices;index_new],[query_subset_indices;index_new])))-(subset_size*max_dist))./(subset_size*subset_size-subset_size);
            if mean_distance_old < mean_distance_new
                break
            end
        end
    end
        

    
if D_old > D_new || (D_old == D_new && mean_distance_old > mean_distance_new)
    index_new = point_old;
end  
    
    %% add new point to subset
    subset_indices = [query_subset_indices;index_new];
    
    %% analysis / plots
    summed_distance = sum(sum(distance_matrix(subset_indices,subset_indices)))-(subset_size*max_dist);
    scatter(ITERATIONS,summed_distance./(subset_size*subset_size-subset_size),'b'), drawnow
    scatter(ITERATIONS, min( min(distance_matrix(subset_indices,subset_indices))),'r'), drawnow
    
end

SUBSET = sort(subset_indices);

%% --- function definitions ---

function [M, I] = permn(V, N, K)
% PERMN - permutations with repetition
%   Using two input variables V and N, M = PERMN(V,N) returns all
%   permutations of N elements taken from the vector V, with repetitions.
%   V can be any type of array (numbers, cells etc.) and M will be of the
%   same type as V.  If V is empty or N is 0, M will be empty.  M has the
%   size numel(V).^N-by-N. 
%
%   When only a subset of these permutations is needed, you can call PERMN
%   with 3 input variables: M = PERMN(V,N,K) returns only the K-ths
%   permutations.  The output is the same as M = PERMN(V,N) ; M = M(K,:),
%   but it avoids memory issues that may occur when there are too many
%   combinations.  This is particulary useful when you only need a few
%   permutations at a given time. If V or K is empty, or N is zero, M will
%   be empty. M has the size numel(K)-by-N. 
%
%   [M, I] = PERMN(...) also returns an index matrix I so that M = V(I).
%
%   Examples:
%     M = permn([1 2 3],2) % returns the 9-by-2 matrix:
%              1     1
%              1     2
%              1     3
%              2     1
%              2     2
%              2     3
%              3     1
%              3     2
%              3     3
%
%     M = permn([99 7],4) % returns the 16-by-4 matrix:
%              99     99    99    99
%              99     99    99     7
%              99     99     7    99
%              99     99     7     7
%              ...
%               7      7     7    99
%               7      7     7     7
%
%     M = permn({'hello!' 1:3},2) % returns the 4-by-2 cell array
%             'hello!'        'hello!'
%             'hello!'        [1x3 double]
%             [1x3 double]    'hello!'
%             [1x3 double]    [1x3 double]
%
%     V = 11:15, N = 3, K = [2 124 21 99]
%     M = permn(V, N, K) % returns the 4-by-3 matrix:
%     %        11  11  12
%     %        15  15  14
%     %        11  15  11
%     %        14  15  14
%     % which are the 2nd, 124th, 21st and 99th permutations
%     % Check with PERMN using two inputs
%     M2 = permn(V,N) ; isequal(M2(K,:),M)
%     % Note that M2 is a 125-by-3 matrix
%
%     % PERMN can be used generate a binary table, as in
%     B = permn([0 1],5)  
%
%   NB Matrix sizes increases exponentially at rate (n^N)*N.
%
%   See also PERMS, NCHOOSEK
%            ALLCOMB, PERMPOS, NEXTPERM, NCHOOSE2 on the File Exchange
% tested in Matlab 2018a
% version 6.2 (jan 2019)
% (c) Jos van der Geest
% Matlab File Exchange Author ID: 10584
% email: samelinoa@gmail.com
% History
% 1.1 updated help text
% 2.0 new faster algorithm
% 3.0 (aug 2006) implemented very fast algorithm
% 3.1 (may 2007) Improved algorithm Roger Stafford pointed out that for some values, the floor
%   operation on floating points, according to the IEEE 754 standard, could return
%   erroneous values. His excellent solution was to add (1/2) to the values
%   of A.
% 3.2 (may 2007) changed help and error messages slightly
% 4.0 (may 2008) again a faster implementation, based on ALLCOMB, suggested on the
%   newsgroup comp.soft-sys.matlab on May 7th 2008 by "Helper". It was
%   pointed out that COMBN(V,N) equals ALLCOMB(V,V,V...) (V repeated N
%   times), ALLCMOB being faster. Actually version 4 is an improvement
%   over version 1 ...
% 4.1 (jan 2010) removed call to FLIPLR, using refered indexing N:-1:1
%   (is faster, suggestion of Jan Simon, jan 2010), removed REPMAT, and
%   let NDGRID handle this
% 4.2 (apr 2011) corrrectly return a column vector for N = 1 (error pointed
%    out by Wilson).
% 4.3 (apr 2013) make a reference to COMBNSUB
% 5.0 (may 2015) NAME CHANGED (COMBN -> PERMN) and updated description,
%   following comment by Stephen Obeldick that this function is misnamed
%   as it produces permutations with repetitions rather then combinations.
% 5.1 (may 2015) always calculate M via indices
% 6.0 (may 2015) merged the functionaly of permnsub (aka combnsub) and this
%   function
% 6.1 (may 2016) fixed spelling errors
% 6.2 (jan 2019) fixed some coding style warnings
narginchk(2, 3) ;
if fix(N) ~= N || N < 0 || numel(N) ~= 1
    error('permn:negativeN','Second argument should be a positive integer') ;
end
nV = numel(V) ;
if nargin==2 
    %% PERMN(V,N) - return all permutations
    if nV == 0 || N == 0
        M = zeros(nV, N) ;
        I = zeros(nV, N) ;
    elseif N == 1
        % return column vectors
        M = V(:) ;
        I = (1:nV).' ;
    else
        % this is faster than the math trick used with 3 inputs below
        [Y{N:-1:1}] = ndgrid(1:nV) ;
        I = reshape(cat(N+1, Y{:}), [], N) ;
        M = V(I) ;
    end
else
    %% PERMN(V,N,K) - return a subset of all permutations
    nK = numel(K) ;
    if nV == 0 || N == 0 || nK == 0
        M = zeros(numel(K), N) ;
        I = zeros(numel(K), N) ;
    elseif nK < 1 || any(K<1) || any(K ~= fix(K))
        error('permn:InvalidIndex','Third argument should contain positive integers.') ;
    else
        V = reshape(V, 1, []) ; % v1.1 make input a row vector
        nV = numel(V) ;
        Npos = nV^N ;
        if any(K > Npos)
            warning('permn:IndexOverflow', ...
                'Values of K exceeding the total number of combinations are saturated.')
            K = min(K, Npos) ;
        end
             
        % The engine is based on version 3.2 with the correction
        % suggested by Roger Stafford. This approach uses a single matrix
        % multiplication.
        B = nV.^(1-N:0) ;
        I = ((K(:)-.5) * B) ; % matrix multiplication
        I = rem(floor(I), nV) + 1 ;
        M = V(I) ;
    end
end
% Algorithm using for-loops
% which can be implemented in C or VB
%
% nv = length(V) ;
% C = zeros(nv^N,N) ; % declaration
% for ii=1:N,
%     cc = 1 ;
%     for jj=1:(nv^(ii-1)),
%         for kk=1:nv,
%             for mm=1:(nv^(N-ii)),
%                 C(cc,ii) = V(kk) ;
%                 cc = cc + 1 ;
%             end
%         end
%     end
% end
end

function [d,A]=strdist(r,b,krk,cas)
%d=strdist(r,b,krk,cas) computes Levenshtein and editor distance 
%between strings r and b with use of Vagner-Fisher algorithm.
%   Levenshtein distance is the minimal quantity of character
%substitutions, deletions and insertions for transformation
%of string r into string b. An editor distance is computed as 
%Levenshtein distance with substitutions weight of 2.
%d=strdist(r) computes numel(r);
%d=strdist(r,b) computes Levenshtein distance between r and b.
%If b is empty string then d=numel(r);
%d=strdist(r,b,krk)computes both Levenshtein and an editor distance
%when krk=2. d=strdist(r,b,krk,cas) computes a distance accordingly 
%with krk and cas. If cas>0 then case is ignored.
%
%Example.
% disp(strdist('matlab'))
%    6
% disp(strdist('matlab','Mathworks'))
%    7
% disp(strdist('matlab','Mathworks',2))
%    7    11
% disp(strdist('matlab','Mathworks',2,1))
%    6     9
switch nargin
   case 1
      d=numel(r);
      return
   case 2
      krk=1;
      bb=b;
      rr=r;
   case 3
       bb=b;
       rr=r;
   case 4
      bb=b;
      rr=r;
      if cas>0
         bb=upper(b);
         rr=upper(r);
      end
end
if krk~=2
   krk=1;
end
d=[];
luma=numel(bb);	lima=numel(rr);
lu1=luma+1;       li1=lima+1;
dl=zeros([lu1,li1]);
dl(1,:)=0:lima;   dl(:,1)=0:luma;
%Distance
for krk1=1:krk
for i=2:lu1
   bbi=bb(i-1);
   for j=2:li1
      kr=krk1;
      if strcmp(rr(j-1),bbi)
         kr=0;
      end
   dl(i,j)=min([dl(i-1,j-1)+kr,dl(i-1,j)+1,dl(i,j-1)+1]);
   end
end
d=[d dl(end,end)];
end


end

function varargout = csvimport( fileName, varargin )
% CSVIMPORT reads the specified CSV file and stores the contents in a cell array or matrix
%
% The file can contain any combination of text & numeric values. Output data format will vary
% depending on the exact composition of the file data.
%
% CSVIMPORT( fileName ):         fileName     -  String specifying the CSV file to be read. Set to
%                                                [] to interactively select the file.
%
% CSVIMPORT( fileName, ... ) : Specify a list of options to be applied when importing the CSV file.
%                              The possible options are:
%                                delimiter     - String to be used as column delimiter. Default
%                                                value is , (comma)
%                                columns       - String or cell array of strings listing the columns
%                                                from which data is to be extracted. If omitted data
%                                                from all columns in the file is imported. If file
%                                                does not contain a header row, the columns
%                                                parameter can be a numeric array listing column
%                                                indices from which data is to be extracted.
%                                outputAsChar  - true / false value indicating whether the data
%                                                should be output as characters. If set to false the
%                                                function attempts to convert each column into a
%                                                numeric array, it outputs the column as characters
%                                                if conversion of any data element in the column
%                                                fails. Default value is false.
%                                uniformOutput - true / false value indicating whether output can be
%                                                returned without encapsulation in a cell array.
%                                                This parameter is ignored if the columns / table
%                                                cannot be converted into a matrix.
%                                noHeader      - true / false value indicating whether the CSV
%                                                file's first line contains column headings. Default
%                                                value is false.
%                                ignoreWSpace  - true / false value indicating whether to ignore
%                                                leading and trailing whitespace in the column
%                                                headers; ignored if noHeader is set to true.
%                                                Default value is false.
%
% The parameters must be specified in the form of param-value pairs, parameter names are not
% case-sensitive and partial matching is supported.
%
% [C1 C2 C3] = CSVIMPORT( fileName, 'columns', {'C1', 'C2', C3'}, ... )
%   This form returns the data from columns in output variables C1, C2 and C3 respectively, the
%   column names are case-sensitive and must match a column name in the file exactly. When fetching
%   data in column mode the number of output columns must match the number of columns to read or it
%   must be one. In the latter case the data from the columns is returned as a single cell matrix.
%
% [C1 C2 C3] = CSVIMPORT( fileName, 'columns', [1, 3, 4], ,'noHeader', true, ... )
%   This form returns the data from columns in output variables C1, C2 and C3 respectively, the
%   columns parameter must contain the column indices when the 'noHeader' option is set to true.
%
% Notes:  1. Function has not been tested on badly formatted CSV files.
%         2. Created using R2007b but has been tested on R2006b.
%
% Revisions:
%   04/28/2009: Corrected typo in an error message
%               Added igonoreWSpace option
%   08/16/2010: Replaced calls to str2num with str2double, the former uses eval leading to unwanted
%               side effects if cells contain text with function names
%
if ( nargin == 0 ) || isempty( fileName )
  [fileName filePath] = uigetfile( '*.csv', 'Select CSV file' );
  if isequal( fileName, 0 )
    return;
  end
  fileName = fullfile( filePath, fileName );
else
  if ~ischar( fileName )
    error( 'csvimport:FileNameError', 'The first argument to %s must be a valid .csv file', ...
      mfilename );
  end
end
%Setup default values
p.delimiter       = ',';
p.columns         = [];
p.outputAsChar    = false;
p.uniformOutput   = true;
p.noHeader        = false;
p.ignoreWSpace    = false;
validParams     = {     ...
  'delimiter',          ...
  'columns',            ...
  'outputAsChar',       ...
  'uniformOutput',      ...
  'noHeader',           ...
  'ignoreWSpace'        ...
  };
%Parse input arguments
if nargin > 1
  if mod( numel( varargin ), 2 ) ~= 0
    error( 'csvimport:InvalidInput', ['All input parameters after the fileName must be in the ' ...
      'form of param-value pairs'] );
  end
  params  = lower( varargin(1:2:end) );
  values  = varargin(2:2:end);
  if ~all( cellfun( @ischar, params ) )
    error( 'csvimport:InvalidInput', ['All input parameters after the fileName must be in the ' ...
      'form of param-value pairs'] );
  end
  lcValidParams   = lower( validParams );
  for ii =  1 : numel( params )
    result        = strmatch( params{ii}, lcValidParams );
    %If unknown param is entered ignore it
    if isempty( result )
      continue
    end
    %If we have multiple matches make sure we don't have a single unambiguous match before throwing
    %an error
    if numel( result ) > 1
      exresult    = strmatch( params{ii}, validParams, 'exact' );
      if ~isempty( exresult )
        result    = exresult;
      else
        %We have multiple possible matches, prompt user to provide an unambiguous match
        error( 'csvimport:InvalidInput', 'Cannot find unambiguous match for parameter ''%s''', ...
          varargin{ii*2-1} );
      end
    end
    result      = validParams{result};
    p.(result)  = values{ii};
  end
end
%Check value attributes
if isempty( p.delimiter ) || ~ischar( p.delimiter )
  error( 'csvimport:InvalidParamType', ['The ''delimiter'' parameter must be a non-empty ' ...
    'character array'] );
end
if isempty( p.noHeader ) || ~islogical( p.noHeader ) || ~isscalar( p.noHeader )
  error( 'csvimport:InvalidParamType', ['The ''noHeader'' parameter must be a non-empty ' ...
    'logical scalar'] );
end
if ~p.noHeader
  if ~isempty( p.columns )
    if ~ischar( p.columns ) && ~iscellstr( p.columns )
      error( 'csvimport:InvalidParamType', ['The ''columns'' parameter must be a character array ' ...
        'or a cell array of strings for CSV files containing column headers on the first line'] );
    end
    if p.ignoreWSpace
      p.columns = strtrim( p.columns );
    end
  end
else
  if ~isempty( p.columns ) && ~isnumeric( p.columns )
    error( 'csvimport:InvalidParamType', ['The ''columns'' parameter must be a numeric array ' ...
      'for CSV files containing column headers on the first line'] );
  end
end
if isempty( p.outputAsChar ) || ~islogical( p.outputAsChar ) || ~isscalar( p.outputAsChar )
  error( 'csvimport:InvalidParamType', ['The ''outputAsChar'' parameter must be a non-empty ' ...
    'logical scalar'] );
end
if isempty( p.uniformOutput ) || ~islogical( p.uniformOutput ) || ~isscalar( p.uniformOutput )
  error( 'csvimport:InvalidParamType', ['The ''uniformOutput'' parameter must be a non-empty ' ...
    'logical scalar'] );
end
%Open file
[fid msg] = fopen( fileName, 'rt' );
if fid == -1
  error( 'csvimport:FileReadError', 'Failed to open ''%s'' for reading.\nError Message: %s', ...
    fileName, msg );
end
colMode         = ~isempty( p.columns );
if ischar( p.columns )
  p.columns     = cellstr( p.columns );
end
nHeaders        = numel( p.columns );
if colMode
  if ( nargout > 1 ) && ( nargout ~= nHeaders )
    error( 'csvimport:NumOutputs', ['The number of output arguments must be 1 or equal to the ' ...
      'number of column names when fetching data for specific columns'] );
  end
end
%Read first line and determine number of columns in data
rowData         = fgetl( fid );
rowData         = regexp( rowData, p.delimiter, 'split' );
nCols           = numel( rowData );
%Check whether all specified columns are present if used in column mode and store their indices
if colMode
  if ~p.noHeader
    if p.ignoreWSpace
      rowData     = strtrim( rowData );
    end
    colIdx        = zeros( 1, nHeaders );
    for ii = 1 : nHeaders
      result      = strmatch( p.columns{ii}, rowData );
      if isempty( result )
        fclose( fid );
        error( 'csvimport:UnknownHeader', ['Cannot locate column header ''%s'' in the file ' ...
          '''%s''. Column header names are case sensitive.'], p.columns{ii}, fileName );
      elseif numel( result ) > 1
        exresult  = strmatch( p.columns{ii}, rowData, 'exact' );
        if numel( exresult ) == 1
          result  = exresult;
        else
          warning( 'csvimport:MultipleHeaderMatches', ['Column header name ''%s'' matched ' ...
            'multiple columns in the file, only the first match (C:%d) will be used.'], ...
            p.columns{ii}, result(1) );
        end
      end
      colIdx(ii)  = result(1);
    end
  else
    colIdx        = p.columns(:);
    if max( colIdx ) > nCols
      fclose( fid );
      error( 'csvimport:BadIndex', ['The specified column index ''%d'' exceeds the number of ' ...
        'columns (%d) in the file'], max( colIdx ), nCols );
    end
  end
end
%Calculate number of lines
pos             = ftell( fid );
if pos == -1
  msg = ferror( fid );
  fclose( fid );
  error( 'csvimport:FileQueryError', 'FTELL on file ''%s'' failed.\nError Message: %s', ...
    fileName, msg );
end
data            = fread( fid );
nLines          = numel( find( data == sprintf( '\n' ) ) ) + 1;
%Reposition file position indicator to beginning of second line
if fseek( fid, pos, 'bof' ) ~= 0
  msg = ferror( fid );
  fclose( fid );
  error( 'csvimport:FileSeekError', 'FSEEK on file ''%s'' failed.\nError Message: %s', ...
    fileName, msg );
end
data            = cell( nLines, nCols );
data(1,:)       = rowData;
emptyRowsIdx    = [];
%Get data for remaining rows
for ii = 2 : nLines
  rowData       = fgetl( fid );
  if isempty( rowData )
    emptyRowsIdx = [emptyRowsIdx(:); ii];
    continue
  end
  rowData       = regexp( rowData, p.delimiter, 'split' );
  nDataElems    = numel( rowData );
  if nDataElems < nCols
    warning( 'csvimport:UnevenColumns', ['Number of data elements on line %d (%d) differs from ' ...
      'that on the first line (%d). Data in this line will be padded.'], ii, nDataElems, nCols );
    rowData(nDataElems+1:nCols) = {''};
  elseif nDataElems > nCols
    warning( 'csvimport:UnevenColumns', ['Number of data elements on line %d (%d) differs from ' ...
      'that one the first line (%d). Data in this line will be truncated.'], ii, nDataElems, nCols );
    rowData     = rowData(1:nCols);
  end
  data(ii,:)    = rowData;
end
%Close file handle
fclose( fid );
data(emptyRowsIdx,:)   = [];
%Process data for final output
uniformOutputPossible  = ~p.outputAsChar;
if p.noHeader
  startRowIdx          = 1;
else
  startRowIdx          = 2;
end
if ~colMode
  if ~p.outputAsChar
    %If we're not outputting the data as characters then try to convert each column to a number
    for ii = 1 : nCols
      colData     = cellfun( @str2double, data(startRowIdx:end,ii), 'UniformOutput', false );
      %If any row contains an entry that cannot be converted to a number then return the whole
      %column as a char array
      if ~any( cellfun( @isnan, colData ) )
        if ~p.noHeader
          data(:,ii)= cat( 1, data(1,ii), colData{:} );
        else
          data(:,ii)= colData;
        end
      end
    end
  end
  varargout{1}    = data;
else
  %In column mode get rid of the headers (if present)
  data            = data(startRowIdx:end,colIdx);
  if ~p.outputAsChar
    %If we're not outputting the data as characters then try to convert each column to a number
    for ii = 1 : nHeaders
      colData     = cellfun( @str2double, data(:,ii), 'UniformOutput', false );
      %If any row contains an entry that cannot be converted to a number then return the whole
      %column as a char array
      if ~any( cellfun( @isnan, colData ) )
        data(:,ii)= colData;
      else
        %If any column cannot be converted to a number then we cannot convert the output to an array
        %or matrix i.e. uniform output is not possible
        uniformOutputPossible = false;
      end
    end
  end
  if nargout == nHeaders
    %Loop through each column and convert to matrix if possible
    for ii = 1 : nHeaders
      if p.uniformOutput && ~any( cellfun( @ischar, data(:,ii) ) )
        varargout{ii} = cell2mat( data(:,ii) );
      else
        varargout{ii} = data(:,ii);
      end
    end
  else
    %Convert entire table to matrix if possible
    if p.uniformOutput && uniformOutputPossible
      data        =  cell2mat( data );
    end
    varargout{1}  = data;
  end
end

end












