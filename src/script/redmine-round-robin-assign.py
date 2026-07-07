#!/usr/bin/python3

# Copyright 2025 IBM, Inc.
# SPDX-License-Identifier: LGPL-2.1-or-later
#
# This script was generated with the assistance of an AI language model.
# AI-Assisted-by: Claude 3.5 Sonnet (Anthropic)
#
# This is free software; you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License version 2.1, as published by
# the Free Software Foundation.  See file COPYING.

"""
Standalone script for round-robin assignment of Redmine issues.

This script assigns unassigned issues to reviewers in a round-robin fashion,
with filtering by project, status, and tags.
"""

import argparse
import logging
import os
import sys
from os.path import expanduser

import redminelib  # https://pypi.org/project/python-redmine/

REDMINE_ENDPOINT = "https://tracker.ceph.com"
REDMINE_API_KEY = None
try:
    with open(expanduser("~/.redmine_key")) as f:
        REDMINE_API_KEY = f.read().strip()
except FileNotFoundError:
    pass
REDMINE_API_KEY = os.getenv("REDMINE_API_KEY", REDMINE_API_KEY)

# Custom field ID for Tags (freeform)
REDMINE_CUSTOM_FIELD_ID_TAGS = 31

log = logging.getLogger(__name__)
log_stream = logging.StreamHandler()
log.addHandler(log_stream)
log.setLevel(logging.INFO)


def round_robin_assign(project, statuses, reviewer_ids, tags=None, dry_run=False, limit=50, debug=False):
    """
    Assign unassigned issues to reviewers in round-robin fashion.
    
    Args:
        project: Redmine project identifier (e.g., "ceph-qa")
        statuses: List of status names to filter (e.g., ["QA Needs Approval"])
        reviewer_ids: List of Redmine user IDs for assignment (numeric)
        tags: Optional list of tags to filter issues
        dry_run: If True, show what would be assigned without making changes
        limit: Maximum number of issues to process
        debug: Enable debug logging
    """
    if debug:
        log.setLevel(logging.DEBUG)
        log.debug("Debug logging enabled.")
    
    if not REDMINE_API_KEY:
        log.fatal("REDMINE_API_KEY not found! Please set REDMINE_API_KEY environment variable or ~/.redmine_key.")
        sys.exit(1)
    
    # Connect to Redmine
    log.info(f"Connecting to {REDMINE_ENDPOINT}")
    R = redminelib.Redmine(REDMINE_ENDPOINT, key=REDMINE_API_KEY)
    log.info("Successfully connected to Redmine.")
    
    # Get project ID
    try:
        log.info(f"Fetching '{project}' project ID from Redmine.")
        proj = R.project.get(project)
        project_id = proj['id']
        log.info(f"Found '{project}' project with ID: {project_id}")
    except redminelib.exceptions.ResourceAttrError:
        log.error(f"Project '{project}' not found in Redmine.")
        sys.exit(1)
    
    # Get status IDs from status names
    log.info("Fetching available statuses from Redmine.")
    all_statuses = R.issue_status.all()
    
    if debug:
        log.debug("Available statuses in Redmine:")
        for status in all_statuses:
            log.debug(f"  - '{status.name}' (ID: {status.id})")
    
    status_ids = []
    for status in all_statuses:
        if status.name in statuses:
            status_ids.append(str(status.id))
            log.info(f"Mapped status '{status.name}' to ID {status.id}")
    
    if not status_ids:
        log.error(f"Could not find any matching status IDs for: {statuses}")
        log.error(f"Available status names: {[s.name for s in all_statuses]}")
        sys.exit(1)
    
    # First, query ALL issues in target statuses to count existing assignments
    log.info("Querying existing assignments in target statuses...")
    all_issues_filter = {
        "project_id": project_id,
        "status_id": ",".join(status_ids),
        "limit": 1000,  # Get more to count assignments accurately
    }
    
    # Add tag filter if specified
    if tags:
        all_issues_filter[f"cf_{REDMINE_CUSTOM_FIELD_ID_TAGS}"] = f"~{tags[0]}"
        log.info(f"Filtering by tag: {tags[0]}")
    
    all_issues = list(R.issue.filter(**all_issues_filter))
    
    # Count current assignments for each reviewer
    current_assignments = {user_id: 0 for user_id in reviewer_ids}
    unassigned_issues = []
    
    for issue in all_issues:
        if hasattr(issue, 'assigned_to') and issue.assigned_to:
            assignee_id = issue.assigned_to.id
            if assignee_id in reviewer_ids:
                current_assignments[assignee_id] += 1
        else:
            unassigned_issues.append(issue)
    
    # Limit unassigned issues to process
    unassigned_issues = unassigned_issues[:limit]
    
    log.info(f"Current assignment counts:")
    for user_id in reviewer_ids:
        log.info(f"  User ID {user_id}: {current_assignments[user_id]} issue(s)")
    
    log.info(f"Found {len(unassigned_issues)} unassigned issue(s) to assign.")
    
    if len(unassigned_issues) == 0:
        log.warning("No unassigned issues found matching the criteria.")
        return
    
    if debug:
        log.debug("Unassigned issues found:")
        for issue in unassigned_issues:
            log.debug(f"  - #{issue.id}: {issue.subject}")
            log.debug(f"    Status: {issue.status.name}")
            log.debug(f"    URL: {REDMINE_ENDPOINT}/issues/{issue.id}")
    
    # Assign issues using load-balanced round-robin
    # Always assign to the reviewer with the fewest current assignments
    new_assignments = {}  # Track new assignments: {user_id: [issue_ids]}
    
    for issue in unassigned_issues:
        # Find reviewer(s) with minimum assignments
        min_count = min(current_assignments.values())
        candidates = [uid for uid in reviewer_ids if current_assignments[uid] == min_count]
        
        # If multiple reviewers have same count, pick the first one in the list
        user_id = candidates[0]
        
        # Update counts
        current_assignments[user_id] += 1
        
        log.debug(f"Processing issue #{issue.id}, assigning to user ID {user_id} (current total: {current_assignments[user_id]})")
        
        try:
            if user_id not in new_assignments:
                new_assignments[user_id] = []
            new_assignments[user_id].append(issue.id)
            
            if dry_run:
                log.info(f"[DRY RUN] Would assign issue #{issue.id} to user ID {user_id}")
            else:
                log.info(f"Assigning issue #{issue.id} to user ID {user_id}")
                try:
                    R.issue.update(issue.id, assigned_to_id=user_id)
                    log.info(f"Successfully assigned issue #{issue.id} to user ID {user_id}")
                except Exception as e:
                    log.error(f"Failed to assign issue #{issue.id}: {e}")
                    if debug:
                        import traceback
                        log.debug(traceback.format_exc())
                    continue
        
        except Exception as e:
            log.error(f"Error processing issue #{issue.id}: {e}")
            if debug:
                import traceback
                log.debug(traceback.format_exc())
            continue
    
    # Print summary
    if new_assignments:
        log.info("\n" + "="*60)
        if dry_run:
            log.info("Load-Balanced Round-Robin Assignment Summary (DRY RUN)")
        else:
            log.info("Load-Balanced Round-Robin Assignment Summary")
        log.info("="*60)
        total_new = sum(len(issues) for issues in new_assignments.values())
        if dry_run:
            log.info(f"Total issues that would be assigned: {total_new}")
        else:
            log.info(f"Total issues assigned: {total_new}")
        log.info("")
        
        # Show final distribution
        log.info("Final assignment distribution:")
        for user_id in reviewer_ids:
            new_count = len(new_assignments.get(user_id, []))
            total_count = current_assignments[user_id]
            if dry_run:
                log.info(f"  User ID {user_id}: {total_count} total ({new_count} new)")
            else:
                log.info(f"  User ID {user_id}: {total_count} total ({new_count} newly assigned)")
        
        log.info("")
        log.info("New assignments:")
        for user_id, issue_ids in sorted(new_assignments.items()):
            log.info(f"  User ID {user_id} - {len(issue_ids)} issue(s):")
            for issue_id in issue_ids:
                log.info(f"    - {REDMINE_ENDPOINT}/issues/{issue_id}")
        log.info("="*60)


def main():
    parser = argparse.ArgumentParser(
        description="Round-robin assignment of Redmine issues",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Dry run to see what would be assigned (use numeric user IDs)
  %(prog)s --project ceph-qa --statuses "QA Needs Approval" \\
    --reviewers 13065,12345 --tags core --dry-run

  # Actually assign issues
  %(prog)s --project ceph-qa --statuses "QA Needs Approval" \\
    --reviewers 13065,12345 --tags core

  # Assign without tag filtering
  %(prog)s --project ceph-qa --statuses "QA Needs Approval" \\
    --reviewers 13065,12345

How to find user IDs:
  1. Go to https://tracker.ceph.com
  2. Click on your username in the top right
  3. Your user ID is in the URL: https://tracker.ceph.com/users/XXXXX
  4. Or search for a user and check their profile URL
        """
    )
    
    parser.add_argument('--project', required=True,
                       help='Redmine project identifier (e.g., "ceph-qa")')
    parser.add_argument('--statuses', required=True,
                       help='Comma-separated list of status names (e.g., "QA Needs Approval,QA Testing")')
    parser.add_argument('--reviewers', required=True,
                       help='Comma-separated list of Redmine user IDs (numeric). '
                            'Find user IDs at https://tracker.ceph.com/users/XXXXX (e.g., "13065,12345")')
    parser.add_argument('--tags',
                       help='Comma-separated list of tags to filter (e.g., "core")')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be assigned without making changes')
    parser.add_argument('--limit', type=int, default=50,
                       help='Maximum number of issues to process (default: 50)')
    parser.add_argument('--debug', action='store_true',
                       help='Enable debug logging')
    
    args = parser.parse_args()
    
    # Parse comma-separated lists
    statuses = [s.strip() for s in args.statuses.split(',') if s.strip()]
    reviewer_ids_str = [r.strip() for r in args.reviewers.split(',') if r.strip()]
    tags = [t.strip() for t in args.tags.split(',') if t.strip()] if args.tags else None
    
    if not statuses:
        log.error("No statuses specified.")
        sys.exit(1)
    
    if not reviewer_ids_str:
        log.error("No reviewer IDs specified.")
        sys.exit(1)
    
    # Convert reviewer IDs to integers
    try:
        reviewer_ids = [int(r) for r in reviewer_ids_str]
    except ValueError as e:
        log.error(f"Invalid reviewer ID format. All reviewer IDs must be numeric.")
        log.error(f"Example: --reviewers '13065,12345'")
        log.error(f"Find user IDs at https://tracker.ceph.com/users/XXXXX")
        sys.exit(1)
    
    round_robin_assign(
        project=args.project,
        statuses=statuses,
        reviewer_ids=reviewer_ids,
        tags=tags,
        dry_run=args.dry_run,
        limit=args.limit,
        debug=args.debug
    )


if __name__ == "__main__":
    main()

# Made with Bob
