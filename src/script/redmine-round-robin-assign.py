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
import random
import re
import sys
from os.path import expanduser
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError
import json
import requests

import redminelib  # https://pypi.org/project/python-redmine/

REDMINE_ENDPOINT = "https://tracker.ceph.com"
REDMINE_API_KEY = None
try:
    with open(expanduser("~/.redmine_key")) as f:
        REDMINE_API_KEY = f.read().strip()
except FileNotFoundError:
    pass
REDMINE_API_KEY = os.getenv("REDMINE_API_KEY", REDMINE_API_KEY)

# GitHub API token (optional, but recommended to avoid rate limits)
GITHUB_TOKEN = None
try:
    with open(expanduser("~/.github_token")) as f:
        GITHUB_TOKEN = f.read().strip()
except FileNotFoundError:
    pass
GITHUB_TOKEN = os.getenv("GITHUB_TOKEN", GITHUB_TOKEN)

# Hardcoded reviewer database
# Maps tracker username -> {tracker_id, github_username, components}
KNOWN_REVIEWERS = {
    "lflores": {
        "tracker_id": 13065,
        "github_username": "ljflores",
        "components": ["core"],
        "slack": "@Laura Flores"
    },
    "nmordech": {
        "tracker_id": 13516,
        "github_username": "NitzanMordhai",
        "components": ["core"],
        "slack": "@Nitzan Mordechai"
    },
    "ksirivad": {
        "tracker_id": 10217,
        "github_username": "kamoltat",
        "components": ["core"],
        "slack": "@Kamoltat (Junior) Sirivadhna"
    },
    "rfriedma": {
        "tracker_id": 10204,
        "github_username": "ronen-fr",
        "components": ["core"],
        "slack": "@Ronen Friedman"
    },
    "amathuri": {
        "tracker_id": 12136,
        "github_username": "amathuria",
        "components": ["core"],
        "slack": "@Aishwarya Mathuria"
    },
    "sseshasa": {
        "tracker_id": 10181,
        "github_username": "sseshasa",
        "components": ["core"],
        "slack": "@Sridhar Seshasayee"
    },
    "shraddhaag": {
        "tracker_id": 11050,
        "github_username": "shraddhaag",
        "components": ["core"],
        "slack": "@Shraddha Agrawal"
    },
    "naveen.naidu@ibm.com": {
        "tracker_id": 15031,
        "github_username": "Naveenaidu",
        "components": ["core"],
        "slack": "@Naveen Naidu"
    },
    "ljsanders": {
        "tracker_id": 14634,
        "github_username": "lee-j-sanders",
        "components": ["core"],
        "slack": "@Lee Sanders"
    },
    "JonBailey1993": {
        "tracker_id": 14872,
        "github_username": "JonBailey1993",
        "components": ["core"],
        "slack": "@Jonathan Bailey"
    },
    "ConnorFawcett": {
        "tracker_id": 14644,
        "github_username": "connorfawcett",
        "components": ["core"],
        "slack": "@Connor Fawcett"
    },
    "Jayaprakash.ceph": {
        "tracker_id": 15114,
        "github_username": "Jayaprakash-ibm",
        "components": ["core"],
        "slack": "@Jaya Prakash"
    },
    # Add more reviewers here as needed
    # "username": {
    #     "tracker_id": 12345,
    #     "github_username": "githubuser",
    #     "components": ["core", "rgw", "rbd"]
    # },
}

def get_reviewer_ids_by_components(components=None):
    """
    Get list of reviewer tracker IDs, optionally filtered by components.
    
    Args:
        components: Optional list of component names to filter reviewers (e.g., ["core", "rgw"])
    
    Returns:
        List of tracker IDs for reviewers matching any of the components
    """
    if components:
        reviewer_ids = []
        for component in components:
            ids = [
                info["tracker_id"]
                for info in KNOWN_REVIEWERS.values()
                if component in info["components"]
            ]
            reviewer_ids.extend(ids)
        # Remove duplicates while preserving order
        seen = set()
        return [x for x in reviewer_ids if not (x in seen or seen.add(x))]
    else:
        return [info["tracker_id"] for info in KNOWN_REVIEWERS.values()]

def get_reviewer_info(username):
    """
    Get reviewer information by tracker username.
    
    Args:
        username: Tracker username (e.g., "lflores")
    
    Returns:
        Dict with tracker_id, github_username, and components, or None if not found
    """
    return KNOWN_REVIEWERS.get(username)

def get_username_by_id(user_id):
    """
    Get tracker username by user ID.
    
    Args:
        user_id: Tracker user ID (numeric)
    
    Returns:
        Username string, or None if not found in known reviewers
    """
    for username, info in KNOWN_REVIEWERS.items():
        if info["tracker_id"] == user_id:
            return username
    return None


def get_slack_id_by_user_id(user_id):
    """
    Get Slack identifier from user ID.
    
    Args:
        user_id: Redmine user ID (numeric)
    
    Returns:
        Slack identifier string if found, None otherwise
    """
    username = get_username_by_id(user_id)
    if username:
        info = get_reviewer_info(username)
        if info and "slack" in info:
            return info["slack"]
    return None


def generate_slack_message(new_assignments, pr_prioritized_assignments, reviewer_ids, dry_run=False):
    """
    Generate a Slack-formatted message for new assignments.
    
    Args:
        new_assignments: Dict mapping user_id -> [issue_ids]
        pr_prioritized_assignments: Dict mapping user_id -> [issue_ids] for PR-prioritized issues
        reviewer_ids: List of all reviewer IDs
        dry_run: Whether this is a dry run
    
    Returns:
        String containing Slack-formatted message
    """
    if not new_assignments:
        return "No new assignments to report."
    
    lines = []
    
    # Header
    if dry_run:
        lines.append("*QA Assignment Summary (DRY RUN)*")
    else:
        lines.append("*QA Assignment Summary*")
    lines.append("")
    
    total_new = sum(len(issues) for issues in new_assignments.values())
    lines.append(f"Total issues assigned: *{total_new}*")
    lines.append("")
    
    # Per-reviewer assignments
    for user_id, issue_ids in sorted(new_assignments.items()):
        slack_id = get_slack_id_by_user_id(user_id)
        username = get_username_by_id(user_id)
        
        if slack_id:
            display_name = slack_id
        elif username:
            display_name = f"@{username}"
        else:
            display_name = f"User ID {user_id}"
        
        pr_count = len(pr_prioritized_assignments.get(user_id, []))
        
        if pr_count > 0:
            lines.append(f"{display_name} - *{len(issue_ids)} issue(s)* ({pr_count} prioritized for PR authorship)")
        else:
            lines.append(f"{display_name} - *{len(issue_ids)} issue(s)*")
        
        for issue_id in issue_ids:
            # Mark PR-prioritized issues with an asterisk
            if user_id in pr_prioritized_assignments and issue_id in pr_prioritized_assignments[user_id]:
                lines.append(f"  • <{REDMINE_ENDPOINT}/issues/{issue_id}|#{issue_id}> *")
            else:
                lines.append(f"  • <{REDMINE_ENDPOINT}/issues/{issue_id}|#{issue_id}>")
        
        lines.append("")
    
    # Add legend if any PR prioritizations occurred
    if pr_prioritized_assignments:
        lines.append("_* = Prioritized due to PR authorship_")
    
    return "\n".join(lines)

def get_reviewer_by_github_username(github_username):
    """
    Get reviewer tracker ID by GitHub username.
    
    Args:
        github_username: GitHub username (e.g., "ljflores")
    
    Returns:
        Tracker user ID, or None if not found in known reviewers
    """
    github_username_lower = github_username.lower()
    for username, info in KNOWN_REVIEWERS.items():
        if info["github_username"].lower() == github_username_lower:
            return info["tracker_id"]
    return None

def get_pr_author_from_github(pr_number):
    """
    Get the author of a GitHub PR using the GitHub API.
    
    Args:
        pr_number: PR number (e.g., 69823)
    
    Returns:
        GitHub username of PR author, or None if unable to fetch
    """
    url = f"https://api.github.com/repos/ceph/ceph/pulls/{pr_number}"
    
    try:
        headers = {'Accept': 'application/vnd.github.v3+json'}
        if GITHUB_TOKEN:
            headers['Authorization'] = f'token {GITHUB_TOKEN}'
        
        req = Request(url, headers=headers)
        with urlopen(req, timeout=5) as response:
            data = json.loads(response.read().decode())
            if 'user' in data and 'login' in data['user']:
                return data['user']['login']
    except (URLError, HTTPError, json.JSONDecodeError, KeyError) as e:
        # Silently fail - we'll just not prioritize this PR
        pass
    
    return None

def check_issue_for_reviewer_prs(issue, reviewer_ids, debug=False):
    """
    Check if any of the PRs in the issue were authored by one of the reviewers.
    
    Args:
        issue: Redmine issue object
        reviewer_ids: List of reviewer user IDs to check
        debug: Enable debug logging
    
    Returns:
        User ID of matching reviewer, or None if no match
    """
    if not hasattr(issue, 'description') or not issue.description:
        return None
    
    # Extract PR numbers from description
    pr_pattern = r'https://github\.com/ceph/ceph/pull/(\d+)'
    pr_numbers = re.findall(pr_pattern, issue.description)
    
    if not pr_numbers:
        return None
    
    if debug:
        log.debug(f"Found {len(pr_numbers)} PR(s) in issue #{issue.id}: {pr_numbers}")
    
    # Build a map of GitHub username -> reviewer user ID
    github_to_reviewer = {}
    for user_id in reviewer_ids:
        username = get_username_by_id(user_id)
        if username:
            info = get_reviewer_info(username)
            if info and info["github_username"]:
                github_to_reviewer[info["github_username"].lower()] = user_id
    
    # Check each PR to see if author matches a reviewer
    for pr_number in pr_numbers:
        pr_author = get_pr_author_from_github(pr_number)
        if pr_author:
            pr_author_lower = pr_author.lower()
            if debug:
                log.debug(f"  PR #{pr_number} author: {pr_author}")
            
            if pr_author_lower in github_to_reviewer:
                matched_user_id = github_to_reviewer[pr_author_lower]
                matched_username = get_username_by_id(matched_user_id)
                if debug:
                    log.debug(f"  Matched reviewer: {matched_username} ({matched_user_id})")
                return matched_user_id
    
    return None

class _IssueWrapper:
    """
    Thin wrapper around a raw Redmine issue JSON dict so the rest of the
    script can use attribute-style access (issue.id, issue.assigned_to.id,
    issue.status.name, issue.description, issue.subject).
    """
    class _Attr:
        def __init__(self, d):
            self.__dict__.update(d)

    def __init__(self, d):
        self.id = d["id"]
        self.subject = d.get("subject", "")
        self.description = d.get("description", "")
        assigned = d.get("assigned_to")
        self.assigned_to = self._Attr(assigned) if assigned else None
        status = d.get("status", {})
        self.status = self._Attr(status)


def fetch_issues(project_id, status_ids, components=None, limit=1000):
    """
    Fetch issues from Redmine using the issue_tags plugin filter, which
    requires the f[]/op[]/v[] envelope that redminelib cannot express for
    plugin-provided fields.

    Args:
        project_id:  Numeric Redmine project ID.
        status_ids:  List of status ID strings.
        components:  Optional list of tag names (e.g. ["core"]).
        limit:       Maximum issues to return.

    Returns:
        List of _IssueWrapper objects.
    """
    url = f"{REDMINE_ENDPOINT}/projects/{project_id}/issues.json"

    # Build the filter parameter list in the order Redmine expects.
    # Repeated keys require a list-of-tuples rather than a dict.
    params = [
        ("key", REDMINE_API_KEY),
        ("f[]", "status_id"),
        ("op[status_id]", "="),
    ]
    for sid in status_ids:
        params.append(("v[status_id][]", sid))

    if components:
        params.append(("f[]", "issue_tags"))
        params.append(("op[issue_tags]", "="))
        for tag in components:
            params.append(("v[issue_tags][]", tag))

    params.append(("f[]", ""))  # close filter list (matches UI behaviour)

    all_issues = []
    offset = 0
    chunk = 100

    while True:
        paged = params + [("limit", chunk), ("offset", offset)]
        resp = requests.get(url, params=paged, timeout=30)
        if resp.status_code != 200:
            raise RuntimeError(
                f"Redmine API returned HTTP {resp.status_code} while fetching issues: {resp.text[:200]}"
            )
        data = resp.json()
        batch = data.get("issues", [])
        all_issues.extend(batch)
        total = data.get("total_count", len(all_issues))
        offset += len(batch)
        if offset >= min(limit, total) or not batch:
            break

    return [_IssueWrapper(d) for d in all_issues[:limit]]


log = logging.getLogger(__name__)
log_stream = logging.StreamHandler()
log.addHandler(log_stream)
log.setLevel(logging.INFO)


def show_distribution_summary(project, statuses, reviewer_ids, components=None, debug=False, slack_format=False):
    """
    Show current distribution of issues among reviewers without assigning.
    
    Args:
        project: Redmine project identifier (e.g., "ceph-qa")
        statuses: List of status names to filter (e.g., ["QA Needs Approval"])
        reviewer_ids: List of Redmine user IDs to check (numeric)
        components: Optional list of components to filter issues
        debug: Enable debug logging
        slack_format: If True, output in Slack-formatted message
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
    
    # Query ALL issues in target statuses
    log.info("Querying issues in target statuses...")
    if components:
        log.info(f"Filtering by component(s): {', '.join(components)}")
    all_issues = fetch_issues(project_id, status_ids, components=components, limit=1000)
    
    # Count current assignments for each reviewer
    current_assignments = {user_id: [] for user_id in reviewer_ids}
    unassigned_issues = []
    other_assignees = {}  # Track issues assigned to non-reviewers
    
    for issue in all_issues:
        if hasattr(issue, 'assigned_to') and issue.assigned_to:
            assignee_id = issue.assigned_to.id
            if assignee_id in reviewer_ids:
                current_assignments[assignee_id].append(issue.id)
            else:
                # Track other assignees
                if assignee_id not in other_assignees:
                    other_assignees[assignee_id] = []
                other_assignees[assignee_id].append(issue.id)
        else:
            unassigned_issues.append(issue)
    
    # Print summary
    if slack_format:
        # Generate Slack-formatted message
        lines = []
        lines.append("*QA Assignment Status*")
        lines.append("")
        
        total_tracked = sum(len(issues) for issues in current_assignments.values())
        lines.append(f"Total issues assigned: *{total_tracked}*")
        if unassigned_issues:
            lines.append(f"Unassigned issues: *{len(unassigned_issues)}*")
        lines.append("")
        
        # Per-reviewer assignments
        for user_id in reviewer_ids:
            issue_ids = current_assignments[user_id]
            if not issue_ids:
                continue  # Skip reviewers with no assignments
            
            slack_id = get_slack_id_by_user_id(user_id)
            username = get_username_by_id(user_id)
            
            if slack_id:
                display_name = slack_id
            elif username:
                display_name = f"@{username}"
            else:
                display_name = f"User ID {user_id}"
            
            lines.append(f"{display_name} - *{len(issue_ids)} issue(s)*")
            for issue_id in issue_ids:
                lines.append(f"  • {REDMINE_ENDPOINT}/issues/{issue_id}")
            lines.append("")
        
        print("\n" + "="*60)
        print("SLACK MESSAGE (copy and paste into Slack):")
        print("="*60)
        print("\n".join(lines))
        print("="*60)
    else:
        # Regular summary output
        log.info("\n" + "="*60)
        log.info("Issue Distribution Summary")
        log.info("="*60)
        log.info(f"Project: {project}")
        log.info(f"Status(es): {', '.join(statuses)}")
        if components:
            log.info(f"Component(s): {', '.join(components)}")
        log.info("")
        
        total_tracked = sum(len(issues) for issues in current_assignments.values())
        log.info(f"Total issues assigned to tracked reviewers: {total_tracked}")
        log.info(f"Total unassigned issues: {len(unassigned_issues)}")
        
        if other_assignees:
            total_other = sum(len(issues) for issues in other_assignees.values())
            log.info(f"Total issues assigned to other users: {total_other}")
        
        log.info("")
        log.info("Distribution among tracked reviewers:")
        for user_id in reviewer_ids:
            count = len(current_assignments[user_id])
            username = get_username_by_id(user_id)
            display_name = f"{username} ({user_id})" if username else f"User ID {user_id}"
            log.info(f"  {display_name}: {count} issue(s)")
            
            if debug and count > 0:
                log.debug(f"    Issues:")
                for issue_id in current_assignments[user_id]:
                    log.debug(f"      - {REDMINE_ENDPOINT}/issues/{issue_id}")
        
        if unassigned_issues:
            log.info("")
            log.info(f"Unassigned issues ({len(unassigned_issues)}):")
            for issue in unassigned_issues[:10]:  # Show first 10
                log.info(f"  - #{issue.id}: {issue.subject}")
                log.info(f"    {REDMINE_ENDPOINT}/issues/{issue.id}")
            if len(unassigned_issues) > 10:
                log.info(f"  ... and {len(unassigned_issues) - 10} more")
        
        log.info("="*60)


def round_robin_assign(project, statuses, reviewer_ids, components=None, dry_run=False, limit=50, debug=False):
    """
    Assign unassigned issues to reviewers in round-robin fashion.
    
    Args:
        project: Redmine project identifier (e.g., "ceph-qa")
        statuses: List of status names to filter (e.g., ["QA Needs Approval"])
        reviewer_ids: List of Redmine user IDs for assignment (numeric)
        components: Optional list of components to filter issues
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
    if components:
        log.info(f"Filtering by component(s): {', '.join(components)}")
    all_issues = fetch_issues(project_id, status_ids, components=components, limit=1000)
    
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
        # Get username for display
        username = get_username_by_id(user_id)
        display_name = f"{username} ({user_id})" if username else f"User ID {user_id}"
        log.info(f"{display_name}: {current_assignments[user_id]} issue(s)")
    
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
    pr_prioritized_assignments = {}  # Track PR-prioritized assignments: {user_id: [issue_ids]}

    # When all reviewers are tied, pick a random starting index so the same
    # person doesn't always get the first issue every run.
    all_tied = len(set(current_assignments.values())) == 1
    if all_tied and len(reviewer_ids) > 1:
        start_index = random.randrange(len(reviewer_ids))
        last_assigned_index = start_index - 1  # rotation picks the next index after this
        log.debug(f"All reviewers tied; randomized start index: {start_index} "
                  f"({get_username_by_id(reviewer_ids[start_index]) or reviewer_ids[start_index]})")
    else:
        last_assigned_index = -1  # default: start from front
    
    for issue in unassigned_issues:
        # Check if any reviewer authored a PR in this issue
        pr_author_id = check_issue_for_reviewer_prs(issue, reviewer_ids, debug=debug)
        
        # Find reviewer(s) with minimum assignments
        min_count = min(current_assignments.values())
        candidates = [uid for uid in reviewer_ids if current_assignments[uid] == min_count]
        
        # Track if this assignment was PR-prioritized
        was_pr_prioritized = False
        
        # If a PR author is among the candidates with minimum count, prioritize them
        if pr_author_id and pr_author_id in candidates:
            user_id = pr_author_id
            last_assigned_index = reviewer_ids.index(user_id)
            username = get_username_by_id(user_id)
            log.info(f"Prioritizing {username} ({user_id}) for issue #{issue.id} (authored PR in this batch)")
            was_pr_prioritized = True
        # If multiple reviewers have same count, rotate through them
        elif len(candidates) > 1:
            # Find where we left off in the reviewer_ids list
            candidate_indices = [reviewer_ids.index(uid) for uid in candidates]
            # Get the next index after last_assigned_index
            next_indices = [idx for idx in candidate_indices if idx > last_assigned_index]
            if next_indices:
                chosen_index = min(next_indices)
            else:
                # Wrap around to the beginning
                chosen_index = min(candidate_indices)
            user_id = reviewer_ids[chosen_index]
            last_assigned_index = chosen_index
        else:
            user_id = candidates[0]
            last_assigned_index = reviewer_ids.index(user_id)
        
        # Update counts
        current_assignments[user_id] += 1
        
        log.debug(f"Processing issue #{issue.id}, assigning to user ID {user_id} (current total: {current_assignments[user_id]})")
        
        try:
            if user_id not in new_assignments:
                new_assignments[user_id] = []
            new_assignments[user_id].append(issue.id)
            
            # Track PR-prioritized assignments separately
            if was_pr_prioritized:
                if user_id not in pr_prioritized_assignments:
                    pr_prioritized_assignments[user_id] = []
                pr_prioritized_assignments[user_id].append(issue.id)
            
            # Get username for display
            username = get_username_by_id(user_id)
            display_name = f"{username} ({user_id})" if username else f"User ID {user_id}"
            
            if dry_run:
                log.info(f"[DRY RUN] Would assign issue #{issue.id} to {display_name}")
            else:
                log.info(f"Assigning issue #{issue.id} to {display_name}")
                try:
                    R.issue.update(issue.id, assigned_to_id=user_id)
                    log.info(f"Successfully assigned issue #{issue.id} to {display_name}")
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
            username = get_username_by_id(user_id)
            display_name = f"{username} ({user_id})" if username else f"User ID {user_id}"
            
            if dry_run:
                log.info(f"  {display_name}: {total_count} total ({new_count} new)")
            else:
                log.info(f"  {display_name}: {total_count} total ({new_count} newly assigned)")
        
        log.info("")
        log.info("New assignments:")
        for user_id, issue_ids in sorted(new_assignments.items()):
            username = get_username_by_id(user_id)
            display_name = f"{username} ({user_id})" if username else f"User ID {user_id}"
            pr_count = len(pr_prioritized_assignments.get(user_id, []))
            
            if pr_count > 0:
                log.info(f"  {display_name} - {len(issue_ids)} issue(s) ({pr_count} prioritized for PR authorship):")
            else:
                log.info(f"  {display_name} - {len(issue_ids)} issue(s):")
            
            for issue_id in issue_ids:
                # Mark PR-prioritized issues with an asterisk
                if user_id in pr_prioritized_assignments and issue_id in pr_prioritized_assignments[user_id]:
                    log.info(f"    - {REDMINE_ENDPOINT}/issues/{issue_id} *")
                else:
                    log.info(f"    - {REDMINE_ENDPOINT}/issues/{issue_id}")
        
        # Add legend if any PR prioritizations occurred
        if pr_prioritized_assignments:
            log.info("")
            log.info("* = Prioritized due to PR authorship")
        
        log.info("="*60)


def main():
    parser = argparse.ArgumentParser(
        description="Round-robin assignment of Redmine issues",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Show current distribution summary (no assignment)
  %(prog)s --project ceph-qa --statuses "QA Needs Approval" \\
    --components core --summary

  # Show current distribution as Slack message (no assignment)
  %(prog)s --project ceph-qa --statuses "QA Needs Approval" \\
    --components core --slack-message

  # Use known reviewers for a component
  %(prog)s --project ceph-qa --statuses "QA Needs Approval" \\
    --components core --dry-run

  # Use known reviewers for multiple components
  %(prog)s --project ceph-qa --statuses "QA Needs Approval" \\
    --components "core,rgw" --dry-run

  # Use specific user IDs
  %(prog)s --project ceph-qa --statuses "QA Needs Approval" \\
    --reviewers 13065,12345 --components core --dry-run

  # Combine known reviewers with additional IDs
  %(prog)s --project ceph-qa --statuses "QA Needs Approval" \\
    --components core --reviewers 99999

  # Assign without component filtering (all unassigned issues)
  %(prog)s --project ceph-qa --statuses "QA Needs Approval" \\
    --reviewers 13065,13516

Known reviewers:
  lflores (13065) - ljflores on GitHub - components: core
  nmordech (13516) - NitzanMordhai on GitHub - components: core

Component matching:
  Issues with multiple components (e.g., "core,rgw") are matched by their
  FIRST component. If --components="core,rgw", an issue tagged "core,rgw"
  will be assigned to a core reviewer.

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
    parser.add_argument('--components',
                       help='Comma-separated list of components to filter and assign (e.g., "core,rgw"). '
                            'Uses known reviewers for these components. Can be combined with --reviewers.')
    parser.add_argument('--reviewers',
                       help='Comma-separated list of Redmine user IDs (numeric). '
                            'Find user IDs at https://tracker.ceph.com/users/XXXXX (e.g., "13065,12345"). '
                            'Can be combined with --components.')
    parser.add_argument('--summary', action='store_true',
                       help='Show current distribution summary without assigning issues. '
                            'Use this flag alone to see how issues are currently distributed.')
    parser.add_argument('--slack-message', action='store_true',
                       help='Show current distribution in Slack-formatted message (uses Slack IDs from known reviewers). '
                            'Like --summary, this queries current state without making assignments.')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be assigned without making changes')
    parser.add_argument('--limit', type=int, default=50,
                       help='Maximum number of issues to process (default: 50)')
    parser.add_argument('--exclude',
                       help='Comma-separated list of tracker usernames to exclude from assignment '
                            '(e.g., "lflores,nmordech")')
    parser.add_argument('--debug', action='store_true',
                       help='Enable debug logging')
    
    args = parser.parse_args()
    
    # Parse comma-separated lists
    statuses = [s.strip() for s in args.statuses.split(',') if s.strip()]
    components = [c.strip() for c in args.components.split(',') if c.strip()] if args.components else None
    
    if not statuses:
        log.error("No statuses specified.")
        sys.exit(1)
    
    # Build reviewer ID list
    reviewer_ids = []
    
    # Add reviewers from components if specified
    if components:
        component_reviewers = get_reviewer_ids_by_components(components)
        if component_reviewers:
            reviewer_ids.extend(component_reviewers)
            log.info(f"Using {len(component_reviewers)} known reviewer(s) for component(s): {', '.join(components)}")
        else:
            log.warning(f"No known reviewers found for component(s): {', '.join(components)}")
    
    # Add additional reviewers from --reviewers if specified
    if args.reviewers:
        reviewer_ids_str = [r.strip() for r in args.reviewers.split(',') if r.strip()]
        try:
            additional_ids = [int(r) for r in reviewer_ids_str]
            reviewer_ids.extend(additional_ids)
            log.info(f"Added {len(additional_ids)} additional reviewer(s) from --reviewers")
        except ValueError as e:
            log.error(f"Invalid reviewer ID format. All reviewer IDs must be numeric.")
            log.error(f"Example: --reviewers '13065,12345'")
            log.error(f"Find user IDs at https://tracker.ceph.com/users/XXXXX")
            sys.exit(1)
    
    # Check that we have at least one reviewer
    if not reviewer_ids:
        log.error("No reviewers specified.")
        log.error("Use --components to select known reviewers, or --reviewers to specify user IDs.")
        log.error("Example: --components core")
        log.error("Example: --components 'core,rgw'")
        log.error("Example: --reviewers '13065,12345'")
        sys.exit(1)
    
    # Remove duplicates while preserving order
    seen = set()
    reviewer_ids = [x for x in reviewer_ids if not (x in seen or seen.add(x))]

    # Exclude specified usernames
    if args.exclude:
        excluded_usernames = [u.strip() for u in args.exclude.split(',') if u.strip()]
        for username in excluded_usernames:
            info = get_reviewer_info(username)
            if info is None:
                log.warning(f"--exclude: unknown username '{username}', skipping.")
                continue
            excluded_id = info["tracker_id"]
            if excluded_id in reviewer_ids:
                reviewer_ids.remove(excluded_id)
                log.info(f"Excluded reviewer '{username}' ({excluded_id}) from assignment.")
            else:
                log.warning(f"--exclude: '{username}' is not in the active reviewer list, skipping.")
    
    # Handle summary/slack-message mode vs assignment mode
    if args.summary or args.slack_message:
        # Summary mode: just show distribution, don't assign
        if args.dry_run:
            log.warning("--dry-run flag is ignored when using --summary or --slack-message")
        if args.limit != 50:
            log.warning("--limit flag is ignored when using --summary or --slack-message")
        
        show_distribution_summary(
            project=args.project,
            statuses=statuses,
            reviewer_ids=reviewer_ids,
            components=components,
            debug=args.debug,
            slack_format=args.slack_message
        )
    else:
        # Assignment mode: perform round-robin assignment
        round_robin_assign(
            project=args.project,
            statuses=statuses,
            reviewer_ids=reviewer_ids,
            components=components,
            dry_run=args.dry_run,
            limit=args.limit,
            debug=args.debug
        )


if __name__ == "__main__":
    main()

# Made with Bob
